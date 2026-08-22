// Streaming, single-pass BLER simulator for the model_out adaptive-m
// selector -- the real thing, not a CSV-and-replay proxy. Modeled directly
// on trace_learnt_path_selector's common/mclass_bler_sim.cpp (same
// decode-once-per-branch, O(1)-memory, --target-errors early-stop,
// --resume-file pattern), extended one more step: for each sample, the
// checkpoint trace features feed a LOADED trace_learnt branch scorer (as
// before) to get the sorted probability vector, which then feeds a SECOND
// loaded model -- the adaptive-m SET-LEVEL selector trained by
// train_adaptive_m.py -- to predict how many of the sorted branches to
// keep. The actual candidate is pickBest (lowest final metric) among that
// top-predicted-m prefix, checked directly against the codeword. No
// per-sample CSV is ever written; mean_m and BLER are running aggregates
// only. This replaces the mclass_adaptive_m_dataset_model_out ->
// plot_cross_t.py CSV-write-then-replay pipeline, which was both far
// slower (full per-sample I/O) and a looser proxy (m_required_basin
// presence-in-prefix, not actually correctness of the pickBest choice).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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
    std::string trace_model_path;
    std::string selector_model_path;
    long channel_seed = -2;
    long message_seed = -1;
    std::string resume_file;
    unsigned int resume_every = 20000;
    unsigned int monitor_every = 2000;
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
        else if (s == "--trace-model") a.trace_model_path = need(s);
        else if (s == "--selector-model") a.selector_model_path = need(s);
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--resume-file") a.resume_file = need(s);
        else if (s == "--resume-every") a.resume_every = (unsigned int)std::stoul(need(s));
        else if (s == "--monitor-every") a.monitor_every = (unsigned int)std::stoul(need(s));
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (!a.snr_set) throw std::runtime_error("--snr is required");
    if (a.samples == 0) throw std::runtime_error("--samples is required");
    if (a.trace_model_path.empty()) throw std::runtime_error("--trace-model is required");
    if (a.selector_model_path.empty()) throw std::runtime_error("--selector-model is required");
    return a;
}

// ---- trace_learnt_path_selector's branch scorer (copied verbatim from
// common/mclass_bler_sim.cpp / mclass_adaptive_m_dataset.cpp) ----
enum class Stat { Active, PmMin, PmGap, PmMean, PmMax, LlrMin, LlrMean, LlrMax };

struct FeatureSlot {
    bool is_branch_onehot = false;
    unsigned int branch_pos = 0;
    unsigned int checkpoint = 0;
    Stat stat = Stat::PmMin;
};

struct TraceModel {
    unsigned int input_dim = 0, hidden = 0, M = 0;
    bool branch_onehot = false, direction_largest = true;
    std::vector<double> mean, std_, W1, b1, W2, b2, W3, b3;
    std::vector<FeatureSlot> slots;
    std::vector<unsigned int> checkpoints_needed;
};

static double statValue(const MClassTracePoint &p, Stat s) {
    switch (s) {
        case Stat::Active: return (double)p.active;
        case Stat::PmMin: return p.pm_min;
        case Stat::PmGap: return p.pm_gap;
        case Stat::PmMean: return p.pm_mean;
        case Stat::PmMax: return p.pm_max;
        case Stat::LlrMin: return p.llr_abs_min;
        case Stat::LlrMean: return p.llr_abs_mean;
        case Stat::LlrMax: return p.llr_abs_max;
    }
    return 0.0;
}

static bool parseFeatureName(const std::string &name, FeatureSlot &slot) {
    if (name.rfind("branch_", 0) == 0) {
        slot.is_branch_onehot = true;
        slot.branch_pos = (unsigned int)std::stoul(name.substr(7));
        return true;
    }
    if (name[0] != 't') return false;
    size_t us = name.find('_');
    if (us == std::string::npos) return false;
    slot.checkpoint = (unsigned int)std::stoul(name.substr(1, us - 1));
    std::string stat = name.substr(us + 1);
    if (stat == "active") slot.stat = Stat::Active;
    else if (stat == "pm_min") slot.stat = Stat::PmMin;
    else if (stat == "pm_gap") slot.stat = Stat::PmGap;
    else if (stat == "pm_mean") slot.stat = Stat::PmMean;
    else if (stat == "pm_max") slot.stat = Stat::PmMax;
    else if (stat == "llr_abs_min") slot.stat = Stat::LlrMin;
    else if (stat == "llr_abs_mean") slot.stat = Stat::LlrMean;
    else if (stat == "llr_abs_max") slot.stat = Stat::LlrMax;
    else return false;
    return true;
}

static std::vector<double> readVec(std::ifstream &in, const std::string &expect_name) {
    std::string name; unsigned int n;
    in >> name >> n;
    if (name != expect_name) throw std::runtime_error("weight file: expected " + expect_name + ", got " + name);
    std::vector<double> v(n);
    for (unsigned int i = 0; i < n; i++) in >> v[i];
    return v;
}

static TraceModel loadTraceModel(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open trace model weights: " + path);
    TraceModel m;
    std::string key;
    in >> key >> m.input_dim;
    in >> key >> m.hidden;
    in >> key >> m.M;
    unsigned int onehot_flag, dir_flag;
    in >> key >> onehot_flag; m.branch_onehot = onehot_flag != 0;
    in >> key >> dir_flag; m.direction_largest = dir_flag != 0;
    m.mean = readVec(in, "mean");
    m.std_ = readVec(in, "std");
    m.W1 = readVec(in, "W1");
    m.b1 = readVec(in, "b1");
    m.W2 = readVec(in, "W2");
    m.b2 = readVec(in, "b2");
    m.W3 = readVec(in, "W3");
    m.b3 = readVec(in, "b3");
    unsigned int n_cols;
    in >> key >> n_cols;
    m.slots.resize(n_cols);
    std::set<unsigned int> cps;
    for (unsigned int i = 0; i < n_cols; i++) {
        std::string name; in >> name;
        if (!parseFeatureName(name, m.slots[i]))
            throw std::runtime_error("cannot parse feature name: " + name);
        if (!m.slots[i].is_branch_onehot) cps.insert(m.slots[i].checkpoint);
    }
    m.checkpoints_needed.assign(cps.begin(), cps.end());
    if (m.slots.size() != m.input_dim)
        throw std::runtime_error("feature_cols count does not match input_dim");
    return m;
}

static double scoreForBranch(const TraceModel &m, unsigned int branch_idx,
                              const std::vector<MClassTracePoint> &traces) {
    std::vector<double> x(m.input_dim);
    for (unsigned int i = 0; i < m.input_dim; i++) {
        const FeatureSlot &s = m.slots[i];
        double raw;
        if (s.is_branch_onehot) {
            raw = (s.branch_pos == branch_idx) ? 1.0 : 0.0;
        } else {
            unsigned int ci = (unsigned int)(std::lower_bound(m.checkpoints_needed.begin(),
                                                                m.checkpoints_needed.end(), s.checkpoint)
                                              - m.checkpoints_needed.begin());
            raw = statValue(traces[ci], s.stat);
        }
        double sd = m.std_[i] < 1e-6 ? 1.0 : m.std_[i];
        x[i] = (raw - m.mean[i]) / sd;
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
    double out = m.b3[0];
    for (unsigned int j = 0; j < m.hidden; j++)
        out += m.W3[j] * h2[j];
    return out;
}

// ---- adaptive-m SET-LEVEL selector, matching
// adaptive_path_selector/*/train_adaptive_m.py's export format exactly:
// input_dim, hidden, log_input, y_mean, y_std, x_mean[N], x_std[N],
// W1,b1,W2,b2,W3,b3 (same Linear-ReLU-Linear-ReLU-Linear shape as the
// trace scorer, but hidden size and log_input are independent per file). ----
struct SelectorModel {
    unsigned int input_dim = 0, hidden = 0;
    bool log_input = true;
    double y_mean = 0.0, y_std = 1.0;
    std::vector<double> x_mean, x_std, W1, b1, W2, b2, W3, b3;
};

static SelectorModel loadSelectorModel(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open selector model weights: " + path);
    SelectorModel m;
    std::string key;
    in >> key >> m.input_dim;
    in >> key >> m.hidden;
    unsigned int log_flag;
    in >> key >> log_flag; m.log_input = log_flag != 0;
    in >> key >> m.y_mean;
    in >> key >> m.y_std;
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

// x_raw = sorted per-branch feature vector (dim input_dim); returns
// predicted m, already clamped to [1, M] and rounded up (ceil), matching
// eval_adaptive_m.py's predict()/score() convention exactly.
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
    double out = m.b3[0];
    for (unsigned int j = 0; j < m.hidden; j++)
        out += m.W3[j] * h2[j];
    double pred = out * m.y_std + m.y_mean;
    long used = (long)std::ceil(pred);
    if (used < 1) used = 1;
    if (used > (long)M) used = M;
    return (unsigned int)used;
}

// ---- resume checkpoint (same shape as mclass_bler_sim.cpp's, adapted for
// this tool's single running total instead of a per-k vector) ----
static std::string configSignature(const Args &a) {
    std::ostringstream ss;
    ss << a.ini_path << "|" << a.snr << "|" << a.samples << "|" << a.target_errors
       << "|" << a.trace_model_path << "|" << a.selector_model_path;
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
        TraceModel trace_model = loadTraceModel(args.trace_model_path);
        SelectorModel selector = loadSelectorModel(args.selector_model_path);
        const unsigned int T = trace_model.checkpoints_needed.back();

        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        const unsigned int n = decoder.codewordLength();
        const unsigned int k = decoder.messageLength();
        const unsigned int M = decoder.branchCount();
        if (T > n) throw std::runtime_error("trace model's max checkpoint exceeds codeword length");
        if (selector.input_dim != M)
            throw std::runtime_error("selector model input_dim does not match branch count M");
        const std::vector<unsigned int> &checkpoints = trace_model.checkpoints_needed;
        const bool largest = trace_model.direction_largest;

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

            // rank by the trace model's own criterion (descending if
            // direction_largest, i.e. basin-style probability)
            std::vector<double> scores(M);
            for (unsigned int b = 0; b < M; b++)
                scores[b] = scoreForBranch(trace_model, b, traces[b]);
            std::vector<unsigned int> order(M);
            for (unsigned int b = 0; b < M; b++) order[b] = b;
            std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
                return largest ? scores[a] > scores[b] : scores[a] < scores[b];
            });

            // selector input: same sorted-probability vector the offline
            // pipeline used (t{T}_trace_prob_sorted_1..M)
            std::vector<double> x(M);
            for (unsigned int r = 0; r < M; r++) {
                double sc = scores[order[r]];
                x[r] = largest ? 1.0 / (1.0 + std::exp(-sc)) : sc;
            }
            unsigned int used_m = selectorPredictM(selector, x, M);
            sum_used_m += used_m;

            // the actual decision: pickBest (lowest final metric) among the
            // top-used_m prefix, checked directly against the codeword --
            // not a m_required_basin proxy.
            unsigned int best = order[0];
            for (unsigned int r = 1; r < used_m; r++)
                if (metrics[order[r]] < metrics[best]) best = order[r];
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
