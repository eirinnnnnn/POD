// Streaming, single-pass Monte Carlo BLER simulator -- modeled directly on
// project/AED/main.cpp's loop (generate -> encode -> noise -> decode ->
// count, O(1) memory, adaptive stop, never materializes per-sample rows to
// disk). Unlike the mclass_dataset -> mclass_trace_dataset -> predict ->
// mclass_eval pipeline, this decodes every branch exactly ONCE per sample
// (decodeBranchTrace already returns the final metric/candidate AND the
// checkpoint trace stats together) and computes full_ped, static-pm-rank
// (k out of M by raw pm_min at checkpoint t), and trace-learned (k out of M
// by a reimplemented MLP forward pass) all from that one decode.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
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
    unsigned int samples = 0;
    unsigned int target_errors = 0;
    std::string stop_metric = "full_ped";  // full_ped | static_trace_min
    unsigned int monitor_every = 0;        // 0 = off; else heartbeat every N samples
    long channel_seed = -2;
    long message_seed = -1;
    unsigned int checkpoint = 0;   // 0 = static/trace disabled, full_ped only
    std::vector<unsigned int> k_list = {8};  // evaluated together from one ranking per sample
    unsigned int stop_k = 0;      // 0 = use max(k_list); which k's count drives the stop condition
    std::string model_path;       // optional -- enables trace-learned
    std::string resume_file;      // optional -- periodic checkpoint so a kill mid-run loses <=resume_every samples, not everything
    unsigned int resume_every = 20000;
};

static std::vector<unsigned int> parseUIntList(const std::string &s) {
    std::vector<unsigned int> out;
    std::stringstream ss(s);
    std::string part;
    while (std::getline(ss, part, ','))
        if (!part.empty()) out.push_back((unsigned int)std::stoul(part));
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

static Args parseArgs(int argc, char **argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string s = argv[i];
        auto need = [&](const std::string &name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value after " + name);
            return argv[++i];
        };
        if (s == "-ini" || s == "--ini") a.ini_path = need(s);
        else if (s == "--snr") a.snr = std::stod(need(s));
        else if (s == "--samples") a.samples = (unsigned int)std::stoul(need(s));
        else if (s == "--target-errors") a.target_errors = (unsigned int)std::stoul(need(s));
        else if (s == "--stop-metric") a.stop_metric = need(s);
        else if (s == "--monitor-every") a.monitor_every = (unsigned int)std::stoul(need(s));
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--checkpoint") a.checkpoint = (unsigned int)std::stoul(need(s));
        else if (s == "--k") a.k_list = {(unsigned int)std::stoul(need(s))};
        else if (s == "--k-list") a.k_list = parseUIntList(need(s));
        else if (s == "--stop-k") a.stop_k = (unsigned int)std::stoul(need(s));
        else if (s == "--model") a.model_path = need(s);
        else if (s == "--resume-file") a.resume_file = need(s);
        else if (s == "--resume-every") a.resume_every = (unsigned int)std::stoul(need(s));
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (a.snr == 0.0) throw std::runtime_error("--snr is required");
    if (a.samples == 0) throw std::runtime_error("--samples is required");
    return a;
}

// ---- MLP scorer (reimplements the trained branch_ranker forward pass) ----
enum class Stat { Active, PmMin, PmGap, PmMean, PmMax, LlrMin, LlrMean, LlrMax };

struct FeatureSlot {
    bool is_branch_onehot = false;
    unsigned int branch_pos = 0;   // which branch index this onehot slot represents
    unsigned int checkpoint = 0;
    Stat stat = Stat::PmMin;
};

struct Model {
    unsigned int input_dim = 0, hidden = 0, M = 0;
    bool branch_onehot = false, direction_largest = true;
    std::vector<double> mean, std_, W1, b1, W2, b2, W3, b3;
    std::vector<FeatureSlot> slots;
    std::vector<unsigned int> checkpoints_needed;  // sorted unique, derived from slots
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

static Model loadModel(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open model weights: " + path);
    Model m;
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

// traces[c] holds the MClassTracePoint for checkpoint m.checkpoints_needed[c] (same order)
static double scoreForBranch(const Model &m, unsigned int branch_idx,
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

// ---- resume checkpoint (periodic, so a kill mid-run loses at most
// resume_every samples instead of the whole run) ----
static std::string configSignature(const Args &a) {
    std::ostringstream ss;
    ss << a.ini_path << "|" << a.snr << "|" << a.samples << "|" << a.checkpoint << "|" << a.model_path << "|k=";
    for (unsigned int kv : a.k_list) ss << kv << ",";
    return ss.str();
}

struct ResumeState {
    bool loaded = false;
    unsigned int sample_count = 0;
    unsigned int full_ped_errors = 0;
    std::vector<unsigned int> static_errors, trace_errors;
    long msg_seed = 0;
    std::string channel_seed_string;
};

static ResumeState loadResume(const std::string &path, const std::string &sig, unsigned int nk) {
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
    in >> key >> st.full_ped_errors;
    in >> key >> st.msg_seed;
    in >> key >> st.channel_seed_string;
    st.static_errors.resize(nk);
    in >> key;
    for (unsigned int i = 0; i < nk; i++) in >> st.static_errors[i];
    st.trace_errors.resize(nk);
    in >> key;
    for (unsigned int i = 0; i < nk; i++) in >> st.trace_errors[i];
    if (!in) {
        std::cerr << "[resume] checkpoint file malformed, starting fresh\n";
        return ResumeState();
    }
    st.loaded = true;
    std::cerr << "[resume] loaded checkpoint at sample " << st.sample_count << "\n";
    return st;
}

static void saveResume(const std::string &path, const std::string &sig, unsigned int sample_count,
                        unsigned int full_ped_errors, long msg_seed, const std::string &channel_seed_string,
                        const std::vector<unsigned int> &static_errors,
                        const std::vector<unsigned int> &trace_errors) {
    if (path.empty()) return;
    std::string tmp = path + ".tmp";
    std::ofstream out(tmp);
    out << "sig " << sig << "\n";
    out << "sample_count " << sample_count << "\n";
    out << "full_ped_errors " << full_ped_errors << "\n";
    out << "msg_seed " << msg_seed << "\n";
    out << "channel_seed_string " << channel_seed_string << "\n";
    out << "static_errors";
    for (unsigned int v : static_errors) out << " " << v;
    out << "\ntrace_errors";
    for (unsigned int v : trace_errors) out << " " << v;
    out << "\n";
    out.close();
    std::rename(tmp.c_str(), path.c_str());  // same filesystem -> atomic
}

// ---- config defaults (same as mclass_dataset.cpp) ----
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
        // NOTE: do NOT force use_AED=true here (unlike mclass_dataset.cpp /
        // mclass_trace_dataset.cpp, which always need AED) -- this tool also
        // runs genuine no-AED SCL-only configs (pure_scl), so use_AED must
        // come from the ini as written. Forcing it true here silently ran
        // every "pure_scl" config as AED M=64 + a big SCL list instead of
        // plain SCL, which is why those results looked implausibly good.
        config["AWGN"]["step_type"] = "SNR";
        config["AWGN"]["start"] = std::to_string(args.snr);
        config["AWGN"]["end"] = std::to_string(args.snr);
        config["AWGN"]["step"] = "1.0";
        config["AWGN"]["seed_string"] = std::to_string(args.channel_seed);

        MClassDecoder decoder(config["AdjustPolarDecoder"]);
        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        Model model;
        bool have_model = !args.model_path.empty();
        if (have_model) model = loadModel(args.model_path);

        std::vector<unsigned int> checkpoints;
        if (args.checkpoint > 0) {
            static const unsigned int ladder[] = {8, 16, 32, 64, 128};
            for (unsigned int c : ladder)
                if (c <= args.checkpoint) checkpoints.push_back(c);
            if (checkpoints.empty() || checkpoints.back() != args.checkpoint)
                checkpoints.push_back(args.checkpoint);
        }
        if (have_model)
            for (unsigned int c : model.checkpoints_needed)
                if (std::find(checkpoints.begin(), checkpoints.end(), c) == checkpoints.end())
                    checkpoints.push_back(c);
        std::sort(checkpoints.begin(), checkpoints.end());

        const unsigned int n = decoder.codewordLength();
        const unsigned int k = decoder.messageLength();
        const unsigned int M = decoder.branchCount();

        long msg_seed = args.message_seed;
        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        const unsigned int nk = (unsigned int)args.k_list.size();
        const unsigned int stop_k = args.stop_k > 0 ? args.stop_k : args.k_list.back();
        unsigned int stop_k_idx = (unsigned int)(std::lower_bound(args.k_list.begin(), args.k_list.end(), stop_k)
                                                   - args.k_list.begin());
        if (stop_k_idx >= nk) stop_k_idx = nk - 1;

        const std::string sig = configSignature(args);
        ResumeState resume = loadResume(args.resume_file, sig, nk);

        unsigned int full_ped_errors = resume.loaded ? resume.full_ped_errors : 0;
        std::vector<unsigned int> static_errors = resume.loaded ? resume.static_errors : std::vector<unsigned int>(nk, 0);
        std::vector<unsigned int> trace_errors = resume.loaded ? resume.trace_errors : std::vector<unsigned int>(nk, 0);
        unsigned int sample_count = resume.loaded ? resume.sample_count : 0;
        const unsigned int start_sample = sample_count;
        unsigned int prev_error_total = 0;
        if (resume.loaded) {
            msg_seed = resume.msg_seed;
            channel.setSeedString(resume.channel_seed_string);
        }

        auto printStatus = [&](unsigned int s) {
            double N = (double)(s + 1);
            std::cerr << "[progress] samples=" << (s + 1)
                       << " full_ped=" << full_ped_errors << " (" << (full_ped_errors / N) << ")";
            if (args.checkpoint > 0) {
                for (unsigned int ki = 0; ki < nk; ki++)
                    std::cerr << " static[k=" << args.k_list[ki] << "]=" << static_errors[ki]
                               << " (" << (static_errors[ki] / N) << ")";
                if (have_model)
                    for (unsigned int ki = 0; ki < nk; ki++)
                        std::cerr << " trace[k=" << args.k_list[ki] << "]=" << trace_errors[ki]
                                   << " (" << (trace_errors[ki] / N) << ")";
            }
            std::cerr << "\n";
            std::cerr.flush();
        };

        for (unsigned int s = start_sample; s < args.samples; s++) {
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

            // full_ped: best metric among all M
            {
                unsigned int best = 0;
                for (unsigned int b = 1; b < M; b++)
                    if (metrics[b] < metrics[best]) best = b;
                if (!correct[best]) full_ped_errors++;
            }

            if (args.checkpoint > 0) {
                unsigned int max_k = std::min(args.k_list.back(), M);

                // static-pm-rank: full sort by pm_min at args.checkpoint (ascending = best),
                // then walk the prefix once, checking correctness at each k breakpoint --
                // top-k for increasing k is just a growing prefix of the same order.
                unsigned int cp_pos = (unsigned int)(std::lower_bound(checkpoints.begin(), checkpoints.end(),
                                                                        args.checkpoint) - checkpoints.begin());
                std::vector<unsigned int> order(M);
                for (unsigned int b = 0; b < M; b++) order[b] = b;
                std::partial_sort(order.begin(), order.begin() + max_k, order.end(),
                                   [&](unsigned int a, unsigned int b) {
                                       return traces[a][cp_pos].pm_min < traces[b][cp_pos].pm_min;
                                   });
                unsigned int best = order[0];
                unsigned int ki = 0;
                for (unsigned int i = 0; i < max_k; i++) {
                    if (i > 0 && metrics[order[i]] < metrics[best]) best = order[i];
                    while (ki < nk && args.k_list[ki] == i + 1) {
                        if (!correct[best]) static_errors[ki]++;
                        ki++;
                    }
                }

                if (have_model) {
                    std::vector<double> scores(M);
                    for (unsigned int b = 0; b < M; b++)
                        scores[b] = scoreForBranch(model, b, traces[b]);
                    std::vector<unsigned int> order2(M);
                    for (unsigned int b = 0; b < M; b++) order2[b] = b;
                    bool largest = model.direction_largest;
                    std::partial_sort(order2.begin(), order2.begin() + max_k, order2.end(),
                                       [&](unsigned int a, unsigned int b) {
                                           return largest ? scores[a] > scores[b] : scores[a] < scores[b];
                                       });
                    unsigned int best2 = order2[0];
                    unsigned int ki2 = 0;
                    for (unsigned int i = 0; i < max_k; i++) {
                        if (i > 0 && metrics[order2[i]] < metrics[best2]) best2 = order2[i];
                        while (ki2 < nk && args.k_list[ki2] == i + 1) {
                            if (!correct[best2]) trace_errors[ki2]++;
                            ki2++;
                        }
                    }
                }
            }

            unsigned int stop_count = full_ped_errors;
            if (args.stop_metric == "static_trace_min")
                stop_count = have_model ? std::min(static_errors[stop_k_idx], trace_errors[stop_k_idx])
                                         : static_errors[stop_k_idx];
            else if (args.stop_metric == "static")
                stop_count = static_errors[stop_k_idx];
            else if (args.stop_metric == "trace")
                stop_count = trace_errors[stop_k_idx];

            unsigned int error_total = full_ped_errors + static_errors[stop_k_idx] + trace_errors[stop_k_idx];
            bool new_error = error_total != prev_error_total;
            prev_error_total = error_total;
            bool heartbeat = args.monitor_every > 0 && ((s + 1) % args.monitor_every == 0);
            if (new_error || heartbeat)
                printStatus(s);

            if (!args.resume_file.empty() && (s + 1) % args.resume_every == 0)
                saveResume(args.resume_file, sig, s + 1, full_ped_errors, msg_seed, channel.getSeedString(),
                           static_errors, trace_errors);

            if (args.target_errors > 0 && stop_count >= args.target_errors)
                break;
        }

        double N = (double)sample_count;
        std::cout << "samples," << sample_count << "\n";
        std::cout << "full_ped_errors," << full_ped_errors << "\n";
        std::cout << "full_ped_bler," << (full_ped_errors / N) << "\n";
        if (args.checkpoint > 0) {
            for (unsigned int ki = 0; ki < nk; ki++) {
                std::cout << "static_errors[k=" << args.k_list[ki] << "]," << static_errors[ki] << "\n";
                std::cout << "static_bler[k=" << args.k_list[ki] << "]," << (static_errors[ki] / N) << "\n";
            }
            if (have_model) {
                for (unsigned int ki = 0; ki < nk; ki++) {
                    std::cout << "trace_errors[k=" << args.k_list[ki] << "]," << trace_errors[ki] << "\n";
                    std::cout << "trace_bler[k=" << args.k_list[ki] << "]," << (trace_errors[ki] / N) << "\n";
                }
            }
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
