// Streaming, single-pass BLER simulator for the pm adaptive-m selector --
// the real thing, not a CSV-and-replay proxy. See
// adaptive_path_selector/model_out/ebch_m7t10/mclass_adaptive_m_bler_sim.cpp
// for the full rationale (same file, mirrored for the pm side): ranks
// branches by raw checkpoint pm_min instead of a trace-model score, no
// trace model needed at all. For each sample: decode all M branches
// through checkpoint t, sort by pm_min, feed the sorted vector to the
// loaded adaptive-m SET-LEVEL selector (train_adaptive_m.py's export
// format) to predict how many to keep, pickBest (lowest final metric)
// among that prefix, check directly against the codeword. No per-sample
// CSV; mean_m and BLER are running aggregates only.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AWGN.h"
#include "libMath.h"
#include "libParser.h"
#include "mclass_decoder_wrapper.h"

struct Args {
    std::string ini_path;
    double snr = 0.0;
    bool snr_set = false;
    unsigned int samples = 0;
    unsigned int target_errors = 0;  // 0 = disabled
    unsigned int checkpoint = 0;
    std::string selector_model_path;
    long channel_seed = -2;
    long message_seed = -1;
    std::string resume_file;
    unsigned int resume_every = 20000;
    unsigned int monitor_every = 2000;
    bool oracle = false;
    unsigned int static_m = 0;  // 0 = disabled
};

static Args parseArgs(int argc, char **argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string s = argv[i];
        auto need = [&](const std::string &name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value after " + name);
            return argv[++i];
        };
        if (s == "-ini" || s == "--ini") a.ini_path = need(s);
        else if (s == "--snr") { a.snr = std::stod(need(s)); a.snr_set = true; }
        else if (s == "--samples") a.samples = (unsigned int)std::stoul(need(s));
        else if (s == "--target-errors") a.target_errors = (unsigned int)std::stoul(need(s));
        else if (s == "--checkpoint") a.checkpoint = (unsigned int)std::stoul(need(s));
        else if (s == "--selector-model") a.selector_model_path = need(s);
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--resume-file") a.resume_file = need(s);
        else if (s == "--resume-every") a.resume_every = (unsigned int)std::stoul(need(s));
        else if (s == "--monitor-every") a.monitor_every = (unsigned int)std::stoul(need(s));
        else if (s == "--oracle") a.oracle = true;
        else if (s == "--static-m") a.static_m = (unsigned int)std::stoul(need(s));
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (!a.snr_set) throw std::runtime_error("--snr is required");
    if (a.samples == 0) throw std::runtime_error("--samples is required");
    if (a.checkpoint == 0) throw std::runtime_error("--checkpoint is required");
    if (a.oracle && a.static_m > 0) throw std::runtime_error("--oracle and --static-m are mutually exclusive");
    if (!a.oracle && a.static_m == 0 && a.selector_model_path.empty())
        throw std::runtime_error("--selector-model is required (unless --oracle or --static-m)");
    return a;
}

// ---- adaptive-m SET-LEVEL selector, matching train_adaptive_m.py's
// export format exactly ----
// Two export formats share this loader:
//  - scalar-regression (train_adaptive_m.py): input_dim, hidden, log_input,
//    y_mean, y_std, x_mean[], x_std[], W1,b1,W2,b2,W3,b3 (W3/b3 size 1) --
//    output_dim implicitly 1.
//  - per-prefix-length BCE classifier (m_prefix/train_adaptive_m_prefix.py):
//    input_dim, hidden, output_dim, log_input, x_mean[], x_std[],
//    W1,b1,W2,b2,W3,b3 (W3/b3 size output_dim==M) -- no y_mean/y_std, output
//    is M logits (one per prefix length k=1..M) instead of one scalar.
struct SelectorModel {
    unsigned int input_dim = 0, hidden = 0, output_dim = 1;
    bool log_input = true;
    double y_mean = 0.0, y_std = 1.0;
    std::vector<double> x_mean, x_std, W1, b1, W2, b2, W3, b3;
};

static std::vector<double> readVec(std::ifstream &in, const std::string &expect_name) {
    std::string name; unsigned int n;
    in >> name >> n;
    if (name != expect_name) throw std::runtime_error("weight file: expected " + expect_name + ", got " + name);
    std::vector<double> v(n);
    for (unsigned int i = 0; i < n; i++) in >> v[i];
    return v;
}

static SelectorModel loadSelectorModel(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open selector model weights: " + path);
    SelectorModel m;
    std::string key;
    in >> key >> m.input_dim;
    in >> key >> m.hidden;
    std::string key3;
    in >> key3;
    if (key3 == "output_dim") {
        in >> m.output_dim;
        unsigned int log_flag;
        in >> key >> log_flag; m.log_input = log_flag != 0;
    } else {
        // key3 == "log_input"
        unsigned int log_flag;
        in >> log_flag; m.log_input = log_flag != 0;
        m.output_dim = 1;
        in >> key >> m.y_mean;
        in >> key >> m.y_std;
    }
    m.x_mean = readVec(in, "x_mean");
    m.x_std = readVec(in, "x_std");
    m.W1 = readVec(in, "W1");
    m.b1 = readVec(in, "b1");
    m.W2 = readVec(in, "W2");
    m.b2 = readVec(in, "b2");
    m.W3 = readVec(in, "W3");
    m.b3 = readVec(in, "b3");
    return m;
}

static unsigned int selectorPredictM(const SelectorModel &m, const std::vector<double> &x_raw, unsigned int M) {
    std::vector<double> x(m.input_dim);
    for (unsigned int i = 0; i < m.input_dim; i++) {
        double v = m.log_input ? std::log1p(x_raw[i]) : x_raw[i];
        double sd = m.x_std[i] < 1e-6 ? 1.0 : m.x_std[i];
        x[i] = (v - m.x_mean[i]) / sd;
    }
    std::vector<double> h1(m.hidden);
    for (unsigned int i = 0; i < m.hidden; i++) {
        double acc = m.b1[i];
        for (unsigned int j = 0; j < m.input_dim; j++)
            acc += m.W1[i * m.input_dim + j] * x[j];
        h1[i] = std::max(0.0, acc);
    }
    std::vector<double> h2(m.hidden);
    for (unsigned int i = 0; i < m.hidden; i++) {
        double acc = m.b2[i];
        for (unsigned int j = 0; j < m.hidden; j++)
            acc += m.W2[i * m.hidden + j] * h1[j];
        h2[i] = std::max(0.0, acc);
    }
    if (m.output_dim == 1) {
        double out = m.b3[0];
        for (unsigned int j = 0; j < m.hidden; j++)
            out += m.W3[j] * h2[j];
        double pred = out * m.y_std + m.y_mean;
        long used = (long)std::ceil(pred);
        if (used < 1) used = 1;
        if (used > (long)M) used = M;
        return (unsigned int)used;
    }

    // per-prefix-length classifier: M logits, one per k=1..M, cumulative-max
    // monotonicity fix, then smallest k crossing p=0.5 (sigmoid(logit)>0
    // equivalently), else clamp to M.
    double running_max = -std::numeric_limits<double>::infinity();
    for (unsigned int k = 0; k < m.output_dim; k++) {
        double logit = m.b3[k];
        for (unsigned int j = 0; j < m.hidden; j++)
            logit += m.W3[k * m.hidden + j] * h2[j];
        running_max = std::max(running_max, logit);
        if (running_max > 0.0)  // sigmoid(x) > 0.5  <=>  x > 0
            return k + 1;
    }
    return M;
}

static std::string configSignature(const Args &a) {
    std::ostringstream ss;
    ss << a.ini_path << "|" << a.snr << "|" << a.samples << "|" << a.target_errors
       << "|" << a.checkpoint << "|"
       << (a.oracle ? "ORACLE" : (a.static_m > 0 ? "STATIC" + std::to_string(a.static_m) : a.selector_model_path));
    return ss.str();
}

struct ResumeState {
    bool loaded = false;
    unsigned int sample_count = 0;
    unsigned int error_count = 0;
    double sum_used_m = 0.0;
    long msg_seed = 0;
    std::string channel_seed_string;
};

static ResumeState loadResume(const std::string &path, const std::string &sig) {
    ResumeState st;
    if (path.empty()) return st;
    std::ifstream in(path);
    if (!in) return st;
    std::string line;
    std::getline(in, line);
    if (line.rfind("sig ", 0) != 0 || line.substr(4) != sig) {
        std::cerr << "[resume] no matching checkpoint, starting fresh\n";
        return st;
    }
    std::string key;
    in >> key >> st.sample_count;
    in >> key >> st.error_count;
    in >> key >> st.sum_used_m;
    in >> key >> st.msg_seed;
    in >> key >> st.channel_seed_string;
    if (!in) {
        std::cerr << "[resume] checkpoint file malformed, starting fresh\n";
        return ResumeState();
    }
    st.loaded = true;
    std::cerr << "[resume] loaded checkpoint at sample " << st.sample_count << "\n";
    return st;
}

static void saveResume(const std::string &path, const std::string &sig, unsigned int sample_count,
                        unsigned int error_count, double sum_used_m, long msg_seed,
                        const std::string &channel_seed_string) {
    if (path.empty()) return;
    std::string tmp = path + ".tmp";
    std::ofstream out(tmp);
    out << "sig " << sig << "\n";
    out << "sample_count " << sample_count << "\n";
    out << "error_count " << error_count << "\n";
    out << "sum_used_m " << sum_used_m << "\n";
    out << "msg_seed " << msg_seed << "\n";
    out << "channel_seed_string " << channel_seed_string << "\n";
    out.close();
    std::rename(tmp.c_str(), path.c_str());
}

static std::map<std::string, std::map<std::string, std::string> > defaultConfig() {
    std::map<std::string, std::map<std::string, std::string> > config;
    config["Monte_Carlo"]["iter_max"] = "10000000";
    config["Monte_Carlo"]["iter_min"] = "0";
    config["Monte_Carlo"]["error_max"] = "50";
    config["Monte_Carlo"]["error_min"] = "50";
    config["Monte_Carlo"]["monitor_slot_size"] = "1";
    config["Monte_Carlo"]["seed_string"] = "-1";
    config["AWGN"]["step_type"] = "SNR";
    config["AWGN"]["start"] = "2.0";
    config["AWGN"]["step"] = "0.25";
    config["AWGN"]["end"] = "7.0";
    config["AWGN"]["seed_string"] = "-2";
    config["AdjustPolarDecoder"]["matrix_src"] = "byGmatrix";
    config["AdjustPolarDecoder"]["Hmatrix_path"] = "";
    config["AdjustPolarDecoder"]["Gmatrix_path"] = "";
    config["AdjustPolarDecoder"]["permutation_src"] = "";
    config["AdjustPolarDecoder"]["permutation_random_seed"] = "-1";
    config["AdjustPolarDecoder"]["target_raw_BER"] = "0.01";
    config["AdjustPolarDecoder"]["operationArray"] = "111011111111";
    config["AdjustPolarDecoder"]["list_size"] = "4";
    config["AdjustPolarDecoder"]["bha_value_setting"] = "????????";
    config["AdjustPolarDecoder"]["OnlyInit"] = "false";
    config["AdjustPolarDecoder"]["use_AED"] = "true";
    config["AdjustPolarDecoder"]["automorphism_src"] = "";
    config["AdjustPolarDecoder"]["aed_L"] = "0";
    return config;
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);

        std::map<std::string, std::map<std::string, std::string> > config = defaultConfig();
        parseConfig(args.ini_path, config);
        config["AdjustPolarDecoder"]["OnlyInit"] = "false";
        config["AWGN"]["step_type"] = "SNR";
        config["AWGN"]["start"] = std::to_string(args.snr);
        config["AWGN"]["end"] = std::to_string(args.snr);
        config["AWGN"]["step"] = "1.0";
        config["AWGN"]["seed_string"] = std::to_string(args.channel_seed);

        MClassDecoder decoder(config["AdjustPolarDecoder"]);
        SelectorModel selector;
        if (!args.oracle && args.static_m == 0) {
            selector = loadSelectorModel(args.selector_model_path);
        }

        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        const unsigned int n = decoder.codewordLength();
        const unsigned int k = decoder.messageLength();
        const unsigned int M = decoder.branchCount();
        if (args.checkpoint > n) throw std::runtime_error("--checkpoint exceeds codeword length");
        if (!args.oracle && args.static_m == 0 && selector.input_dim != M)
            throw std::runtime_error("selector model input_dim does not match branch count M");

        std::vector<unsigned int> checkpoints;
        { static const unsigned int ladder[] = {8, 16, 32, 64, 128};
          for (unsigned int c : ladder) if (c <= args.checkpoint) checkpoints.push_back(c);
          if (checkpoints.empty() || checkpoints.back() != args.checkpoint) checkpoints.push_back(args.checkpoint); }
        const unsigned int cp_pos = (unsigned int)(std::lower_bound(checkpoints.begin(), checkpoints.end(),
                                                                      args.checkpoint) - checkpoints.begin());

        const std::string sig = configSignature(args);
        ResumeState resume = loadResume(args.resume_file, sig);

        unsigned int error_count = resume.loaded ? resume.error_count : 0;
        double sum_used_m = resume.loaded ? resume.sum_used_m : 0.0;
        unsigned int sample_count = resume.loaded ? resume.sample_count : 0;
        const unsigned int start_sample = sample_count;
        long msg_seed = resume.loaded ? resume.msg_seed : args.message_seed;
        if (resume.loaded) channel.setSeedString(resume.channel_seed_string);

        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        auto printStatus = [&](unsigned int s) {
            double N = (double)(s + 1);
            std::cerr << "[progress] samples=" << (s + 1) << " errors=" << error_count
                       << " (" << (error_count / N) << ") mean_m=" << (sum_used_m / N) << "\n";
            std::cerr.flush();
        };

        unsigned int s = start_sample;
        for (; s < args.samples; s++) {
            for (unsigned int i = 0; i < k; i++)
                message[i] = (ran0(&msg_seed) > 0.5 ? 1 : 0);
            decoder.doEncode(message, codeword);
            channel.addNoise(codeword, received);

            std::vector<double> metrics(M, 0.0);
            std::vector<char> correct(M, 0);
            std::vector<std::vector<MClassTracePoint> > traces(M);
            for (unsigned int b = 0; b < M; b++) {
                std::vector<char> candidate;
                double metric = 0.0;
                bool ok = decoder.decodeBranchTrace(received, b, checkpoints, candidate, metric, traces[b]);
                metrics[b] = ok ? metric : std::numeric_limits<double>::max();
                correct[b] = ok && (candidate == codeword);
            }
            sample_count = s + 1;

            std::vector<unsigned int> order(M);
            for (unsigned int b = 0; b < M; b++) order[b] = b;
            std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
                return traces[a][cp_pos].pm_min < traces[b][cp_pos].pm_min;
            });

            unsigned int used_m;
            unsigned int best = order[0];
            if (args.oracle) {
                // true per-sample minimal prefix length under THIS ranking:
                // walk the running-best ratchet and stop at the first r
                // where it's already correct (pickBest's own oracle) -- no
                // selector model involved. Only fails when full_ped (r=M)
                // itself fails, i.e. this is the real floor.
                used_m = M;
                for (unsigned int r = 1; r <= M; r++) {
                    if (r > 1 && metrics[order[r - 1]] < metrics[best]) best = order[r - 1];
                    if (correct[best]) { used_m = r; break; }
                }
            } else if (args.static_m > 0) {
                // fixed, non-adaptive prefix length for every sample --
                // isolates the value of per-sample adaptivity from just
                // picking a good constant width.
                used_m = args.static_m;
                for (unsigned int r = 1; r < used_m; r++)
                    if (metrics[order[r]] < metrics[best]) best = order[r];
            } else {
                std::vector<double> x(M);
                for (unsigned int r = 0; r < M; r++)
                    x[r] = traces[order[r]][cp_pos].pm_min;
                used_m = selectorPredictM(selector, x, M);
                for (unsigned int r = 1; r < used_m; r++)
                    if (metrics[order[r]] < metrics[best]) best = order[r];
            }
            sum_used_m += used_m;
            if (!correct[best]) error_count++;

            if ((s + 1) % args.monitor_every == 0) printStatus(s);
            if (args.resume_file.size() && (s + 1) % args.resume_every == 0)
                saveResume(args.resume_file, sig, sample_count, error_count, sum_used_m,
                           msg_seed, channel.getSeedString());
            if (args.target_errors > 0 && error_count >= args.target_errors) {
                sample_count = s + 1;
                break;
            }
        }

        std::cout << "samples," << sample_count << "\n";
        std::cout << "errors," << error_count << "\n";
        std::cout << "bler," << ((double)error_count / (double)sample_count) << "\n";
        std::cout << "mean_m," << (sum_used_m / (double)sample_count) << "\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
