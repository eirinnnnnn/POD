// Streaming, single-pass Monte Carlo BLER simulator for the Golay(24,12)
// code, modeled on project/MClassifier/common/mclass_bler_sim.cpp.
//
// Genuine early-path-selection: branches are ranked for pruning by a
// CHECKPOINT (partial, mid-decode) trace snapshot at decode_idx=t (out of
// codeword_length=24), then the top-m survivors are run to completion and
// the best among them (prefer parity-valid, lowest FINAL metric) is
// selected. Because the ranking criterion (checkpoint signal) and the
// selection criterion (final metric) are different quantities, this can
// genuinely diverge from full-PED -- unlike ranking and selecting by the
// same final metric, which provably cannot (see PolarMultiKernal_PS_PED's
// decodeAllBranchesTrace()/decodeOnceWithOrderTrace()). Two ranking
// strategies are supported per (checkpoint) snapshot: static-pm-rank
// (raw checkpoint pm_min, no model) and trace-learned (an MLP over the
// checkpoint's PSPEDTracePoint stats + branch one-hot, trained by
// architectures/trace/train_mclass_trace_ranker.py on golay_trace_dataset
// output and exported by architectures/trace/export_model_weights.py --
// both reused unchanged from the eBCH pipeline since the weight file
// format and feature-name convention are code-agnostic).
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
#include "PolarMultiKernal_PS_PED.h"

struct Args {
    std::string ini_path;
    double snr = 0.0;
    bool snr_set = false;
    unsigned int samples = 0;
    unsigned int target_errors = 0;
    unsigned int monitor_every = 0;
    long channel_seed = -2;
    long message_seed = -1;
    unsigned int checkpoint = 0;              // decode_idx progress count (1..codeword_length) used to rank branches
    std::vector<unsigned int> m_list = {8};   // evaluated together from one ranking per sample
    unsigned int stop_m = 0;                  // 0 = use max(m_list)
    std::string model_path;                   // optional -- enables trace-learned
    std::string resume_file;
    unsigned int resume_every = 20000;
    bool dump_frozen = false;                  // diagnostic: print diverge_flag (info vs frozen per decode_idx) and exit
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
        else if (s == "--snr") { a.snr = std::stod(need(s)); a.snr_set = true; }
        else if (s == "--samples") a.samples = (unsigned int)std::stoul(need(s));
        else if (s == "--target-errors") a.target_errors = (unsigned int)std::stoul(need(s));
        else if (s == "--monitor-every") a.monitor_every = (unsigned int)std::stoul(need(s));
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--checkpoint") a.checkpoint = (unsigned int)std::stoul(need(s));
        else if (s == "--m") a.m_list = {(unsigned int)std::stoul(need(s))};
        else if (s == "--m-list") a.m_list = parseUIntList(need(s));
        else if (s == "--stop-m") a.stop_m = (unsigned int)std::stoul(need(s));
        else if (s == "--model") a.model_path = need(s);
        else if (s == "--resume-file") a.resume_file = need(s);
        else if (s == "--resume-every") a.resume_every = (unsigned int)std::stoul(need(s));
        else if (s == "--dump-frozen") a.dump_frozen = true;
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (!a.snr_set) throw std::runtime_error("--snr is required");
    if (a.samples == 0) throw std::runtime_error("--samples is required");
    if (a.checkpoint == 0) throw std::runtime_error("--checkpoint is required (t, decode_idx progress count)");
    return a;
}

// ---- MLP scorer (reimplements the trained branch_ranker forward pass) ----
enum class Stat { Active, PmMin, PmGap, PmMean, PmMax, LlrMin, LlrMean, LlrMax };

struct FeatureSlot {
    bool is_branch_onehot = false;
    unsigned int branch_pos = 0;
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

static double statValue(const PSPEDTracePoint &p, Stat s) {
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

// traces[c] holds the PSPEDTracePoint for checkpoint m.checkpoints_needed[c] (same order)
static double scoreForBranch(const Model &m, unsigned int branch_idx,
                              const std::vector<PSPEDTracePoint> &traces) {
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

// ---- resume checkpoint ----
static std::string configSignature(const Args &a) {
    std::ostringstream ss;
    ss << a.ini_path << "|" << a.snr << "|" << a.samples << "|t=" << a.checkpoint << "|" << a.model_path << "|m=";
    for (unsigned int mv : a.m_list) ss << mv << ",";
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

static ResumeState loadResume(const std::string &path, const std::string &sig, unsigned int nm) {
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
    st.static_errors.resize(nm);
    in >> key;
    for (unsigned int i = 0; i < nm; i++) in >> st.static_errors[i];
    st.trace_errors.resize(nm);
    in >> key;
    for (unsigned int i = 0; i < nm; i++) in >> st.trace_errors[i];
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
    std::rename(tmp.c_str(), path.c_str());
}

// ---- config defaults (mirrors main_AED.cpp) ----
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
    config["PolarMultiKernal"]["matrix_src"] = "byGmatrix";
    config["PolarMultiKernal"]["Hmatrix_path"] = "";
    config["PolarMultiKernal"]["Gmatrix_path"] = "";
    config["PolarMultiKernal"]["permutation_src"] = "random";
    config["PolarMultiKernal"]["permutation_random_seed"] = "-1";
    config["PolarMultiKernal"]["dynamic_frozen_process"] = "frozen";
    config["PolarMultiKernal"]["target_raw_BER"] = "0.01";
    config["PolarMultiKernal"]["list_size"] = "4";
    config["PolarMultiKernal"]["kernal_string"] = "23,23,23,753";
    config["PolarMultiKernal"]["bha_value_setting"] = "";
    config["PolarMultiKernal"]["message_length"] = "12";
    config["PolarMultiKernal"]["use_AED"] = "true";
    config["PolarMultiKernal"]["automorphism_src"] = "";
    config["PolarMultiKernal"]["aed_L"] = "0";
    return config;
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);

        std::map<std::string, std::map<std::string, std::string> > config = defaultConfig();
        parseConfig(args.ini_path, config);
        config["AWGN"]["step_type"] = "SNR";
        config["AWGN"]["start"] = std::to_string(args.snr);
        config["AWGN"]["end"] = std::to_string(args.snr);
        config["AWGN"]["step"] = "1.0";
        config["AWGN"]["seed_string"] = std::to_string(args.channel_seed);

        PolarMultiKernal decoder(config["PolarMultiKernal"]);

        if (args.dump_frozen) {
            const std::vector<char> &df = decoder.divergeFlags();
            for (unsigned int i = 0; i < df.size(); i++)
                std::cout << "decode_idx=" << i << " " << (df[i] ? "info" : "frozen")
                          << " relation_span=" << decoder.relationSpan(i) << "\n";
            return 0;
        }

        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        Model model;
        bool have_model = !args.model_path.empty();
        if (have_model) model = loadModel(args.model_path);

        const unsigned int n = decoder.getCodewordLength();
        const unsigned int k = decoder.getMessageLength();
        if (args.checkpoint > n) throw std::runtime_error("--checkpoint exceeds codeword length");

        // Cumulative checkpoint ladder up to args.checkpoint (matches the
        // eBCH "upto_tX" convention), plus any additional checkpoints the
        // trace-learned model needs beyond that ladder.
        std::vector<unsigned int> checkpoints;
        static const unsigned int ladder[] = {4, 8, 16, 20};
        for (unsigned int c : ladder)
            if (c <= args.checkpoint) checkpoints.push_back(c);
        if (checkpoints.empty() || checkpoints.back() != args.checkpoint)
            checkpoints.push_back(args.checkpoint);
        if (have_model)
            for (unsigned int c : model.checkpoints_needed)
                if (std::find(checkpoints.begin(), checkpoints.end(), c) == checkpoints.end())
                    checkpoints.push_back(c);
        std::sort(checkpoints.begin(), checkpoints.end());
        const unsigned int rank_cp_pos = (unsigned int)(std::lower_bound(checkpoints.begin(), checkpoints.end(),
                                                                          args.checkpoint) - checkpoints.begin());

        long msg_seed = args.message_seed;
        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        const unsigned int nm = (unsigned int)args.m_list.size();
        const unsigned int stop_m = args.stop_m > 0 ? args.stop_m : args.m_list.back();
        unsigned int stop_m_idx = (unsigned int)(std::lower_bound(args.m_list.begin(), args.m_list.end(), stop_m)
                                                   - args.m_list.begin());
        if (stop_m_idx >= nm) stop_m_idx = nm - 1;

        const std::string sig = configSignature(args);
        ResumeState resume = loadResume(args.resume_file, sig, nm);

        unsigned int full_ped_errors = resume.loaded ? resume.full_ped_errors : 0;
        std::vector<unsigned int> static_errors = resume.loaded ? resume.static_errors : std::vector<unsigned int>(nm, 0);
        std::vector<unsigned int> trace_errors = resume.loaded ? resume.trace_errors : std::vector<unsigned int>(nm, 0);
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
            for (unsigned int mi = 0; mi < nm; mi++)
                std::cerr << " static[m=" << args.m_list[mi] << "]=" << static_errors[mi]
                           << " (" << (static_errors[mi] / N) << ")";
            if (have_model)
                for (unsigned int mi = 0; mi < nm; mi++)
                    std::cerr << " trace[m=" << args.m_list[mi] << "]=" << trace_errors[mi]
                               << " (" << (trace_errors[mi] / N) << ")";
            std::cerr << "\n";
            std::cerr.flush();
        };

        std::string log;
        for (unsigned int s = start_sample; s < args.samples; s++) {
            for (unsigned int i = 0; i < k; i++)
                message[i] = (ran0(&msg_seed) > 0.5 ? 1 : 0);
            decoder.doEncode(message, codeword);
            channel.addNoise(codeword, received);

            std::vector<std::vector<char> > candidates;
            std::vector<double> metrics;
            std::vector<char> valid;
            std::vector<std::vector<PSPEDTracePoint> > checkpoint_trace;  // [branch][0]
            decoder.decodeAllBranchesTrace(received, candidates, metrics, valid, checkpoints, checkpoint_trace, log);
            const unsigned int L = (unsigned int)candidates.size();

            std::vector<char> correct(L, 0);
            for (unsigned int b = 0; b < L; b++)
                correct[b] = (candidates[b] == codeword);

            // prefer a parity-valid candidate with lowest FINAL metric,
            // fall back to lowest final metric overall (mirrors doDecode()'s policy)
            auto pickBest = [&](const std::vector<unsigned int> &idxs) {
                int best_valid = -1, best_any = -1;
                for (unsigned int idx : idxs) {
                    if (best_any < 0 || metrics[idx] < metrics[(unsigned int)best_any]) best_any = (int)idx;
                    if (valid[idx] && (best_valid < 0 || metrics[idx] < metrics[(unsigned int)best_valid]))
                        best_valid = (int)idx;
                }
                return best_valid >= 0 ? (unsigned int)best_valid : (unsigned int)best_any;
            };

            {
                std::vector<unsigned int> all(L);
                for (unsigned int b = 0; b < L; b++) all[b] = b;
                unsigned int best = pickBest(all);
                if (!correct[best]) full_ped_errors++;
            }

            unsigned int max_m = std::min(args.m_list.back(), L);

            // static-pm-rank: sort all L branches by their CHECKPOINT pm_min
            // at args.checkpoint (ascending), keep the top-m survivors,
            // then pick the best among that pruned subset by FINAL metric
            // (prefer valid).
            {
                std::vector<unsigned int> order(L);
                for (unsigned int b = 0; b < L; b++) order[b] = b;
                std::partial_sort(order.begin(), order.begin() + max_m, order.end(),
                                   [&](unsigned int a, unsigned int b) {
                                       return checkpoint_trace[a][rank_cp_pos].pm_min
                                            < checkpoint_trace[b][rank_cp_pos].pm_min;
                                   });
                std::vector<unsigned int> prefix;
                unsigned int mi = 0;
                for (unsigned int i = 0; i < max_m; i++) {
                    prefix.push_back(order[i]);
                    while (mi < nm && args.m_list[mi] == i + 1) {
                        unsigned int best = pickBest(prefix);
                        if (!correct[best]) static_errors[mi]++;
                        mi++;
                    }
                }
            }

            // trace-learned: same idea, but rank by the trained model's
            // score instead of raw pm_min.
            if (have_model) {
                std::vector<double> scores(L);
                for (unsigned int b = 0; b < L; b++)
                    scores[b] = scoreForBranch(model, b, checkpoint_trace[b]);
                std::vector<unsigned int> order2(L);
                for (unsigned int b = 0; b < L; b++) order2[b] = b;
                bool largest = model.direction_largest;
                std::partial_sort(order2.begin(), order2.begin() + max_m, order2.end(),
                                   [&](unsigned int a, unsigned int b) {
                                       return largest ? scores[a] > scores[b] : scores[a] < scores[b];
                                   });
                std::vector<unsigned int> prefix2;
                unsigned int mi2 = 0;
                for (unsigned int i = 0; i < max_m; i++) {
                    prefix2.push_back(order2[i]);
                    while (mi2 < nm && args.m_list[mi2] == i + 1) {
                        unsigned int best2 = pickBest(prefix2);
                        if (!correct[best2]) trace_errors[mi2]++;
                        mi2++;
                    }
                }
            }
            sample_count = s + 1;

            unsigned int stop_count = have_model ? std::min(static_errors[stop_m_idx], trace_errors[stop_m_idx])
                                                  : static_errors[stop_m_idx];
            unsigned int error_total = full_ped_errors + static_errors[stop_m_idx] + trace_errors[stop_m_idx];
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
        for (unsigned int mi = 0; mi < nm; mi++) {
            std::cout << "static_errors[m=" << args.m_list[mi] << "]," << static_errors[mi] << "\n";
            std::cout << "static_bler[m=" << args.m_list[mi] << "]," << (static_errors[mi] / N) << "\n";
        }
        if (have_model) {
            for (unsigned int mi = 0; mi < nm; mi++) {
                std::cout << "trace_errors[m=" << args.m_list[mi] << "]," << trace_errors[mi] << "\n";
                std::cout << "trace_bler[m=" << args.m_list[mi] << "]," << (trace_errors[mi] / N) << "\n";
            }
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
