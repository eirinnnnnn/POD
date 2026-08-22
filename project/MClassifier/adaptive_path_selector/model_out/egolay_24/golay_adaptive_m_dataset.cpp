// Dataset generator for the "adaptive path selector, model_out variant"
// (Golay24 side): same SET-LEVEL target as adaptive_path_selector/pm
// (predict, per sample, the minimum pruning width m such that keeping the
// top-m branches still contains the branch full-PED itself would have
// picked), but ranked and featurized by the trace_learnt_path_selector's
// OWN per-branch output instead of the raw checkpoint pm_min.
//
// trace_learnt_path_selector's model IO (see
// trace_learnt_path_selector/egolay_24/common/golay_mclass_bler_sim.cpp):
// input per branch is that branch's path-metric AND channel-LLR-magnitude
// statistics (pm_min/pm_gap/pm_mean/pm_max, llr_abs_min/mean/max) at every
// checkpoint up to t, plus a branch-identity one-hot; output is a single
// scalar per branch (a BCEWithLogits logit for target=basin models --
// sigmoid(logit) is then a literal per-branch correctness probability).
// Stacking that scalar over all L branches gives exactly the dim-L
// probability vector this generator uses as its own input feature,
// replacing pm's t{T}_pm_min_sorted_1..L columns with
// t{T}_trace_prob_sorted_1..L (sorted by the trace model's own ranking
// criterion, not by pm_min). T is taken from the trace model itself (its
// highest trained checkpoint), not a separate --checkpoint argument.
//
// This is a NEW file only -- copies (does not modify) the Model/loadModel/
// scoreForBranch block from trace_learnt_path_selector's
// common/golay_mclass_bler_sim.cpp, and keeps the SAME teacher convention
// as adaptive_path_selector/pm/egolay_24/golay_adaptive_m_dataset.cpp
// (pickBest: prefer a parity-valid candidate with the lowest FINAL metric,
// else lowest final metric overall -- matches golay_mclass_bler_sim.cpp's
// own full_ped reference, NOT eBCH's plain-argmin convention).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
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
    std::string trace_model_path;
    long channel_seed = -2;
    long message_seed = -1;
    std::string out_path;
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
        else if (s == "--trace-model") a.trace_model_path = need(s);
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--out") a.out_path = need(s);
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (!a.snr_set) throw std::runtime_error("--snr is required");
    if (a.samples == 0) throw std::runtime_error("--samples is required");
    if (a.trace_model_path.empty()) throw std::runtime_error("--trace-model is required");
    if (a.out_path.empty()) throw std::runtime_error("--out is required");
    return a;
}

// ---- trace_learnt_path_selector's branch scorer, copied verbatim from
// common/golay_mclass_bler_sim.cpp (Model/loadModel/scoreForBranch) ----
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
    std::vector<unsigned int> checkpoints_needed;
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
        Model model = loadModel(args.trace_model_path);
        const unsigned int T = model.checkpoints_needed.back();

        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        const unsigned int n = decoder.getCodewordLength();
        const unsigned int k = decoder.getMessageLength();
        if (T > n) throw std::runtime_error("trace model's max checkpoint exceeds codeword length");
        const std::vector<unsigned int> &checkpoints = model.checkpoints_needed;
        const bool largest = model.direction_largest;
        const char *col_base = largest ? "trace_prob_sorted" : "trace_score_sorted";

        long msg_seed = args.message_seed;
        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        const unsigned int L_fixed = decoder.branchCount();

        std::ofstream out(args.out_path);
        if (!out) throw std::runtime_error("cannot open --out for writing: " + args.out_path);
        out << "sample,L,m_required,full_ped_correct";
        for (unsigned int r = 1; r <= L_fixed; r++)
            out << ",t" << T << "_" << col_base << "_" << r;
        out << "\n";

        std::string log;

        for (unsigned int s = 0; s < args.samples; s++) {
            for (unsigned int i = 0; i < k; i++)
                message[i] = (ran0(&msg_seed) > 0.5 ? 1 : 0);
            decoder.doEncode(message, codeword);
            channel.addNoise(codeword, received);

            std::vector<std::vector<char> > candidates;
            std::vector<double> metrics;
            std::vector<char> valid;
            std::vector<std::vector<PSPEDTracePoint> > checkpoint_trace;
            decoder.decodeAllBranchesTrace(received, candidates, metrics, valid, checkpoints, checkpoint_trace, log);
            const unsigned int L = (unsigned int)candidates.size();
            if (L != L_fixed) throw std::runtime_error("branch count L changed mid-run (unexpected)");

            std::vector<char> correct(L, 0);
            for (unsigned int b = 0; b < L; b++)
                correct[b] = (candidates[b] == codeword);

            // teacher = pickBest(all L branches) -- same convention as
            // adaptive_path_selector/pm/egolay_24 and
            // golay_mclass_bler_sim.cpp's full_ped reference.
            int best_valid = -1, best_any = -1;
            for (unsigned int b = 0; b < L; b++) {
                if (best_any < 0 || metrics[b] < metrics[(unsigned int)best_any]) best_any = (int)b;
                if (valid[b] && (best_valid < 0 || metrics[b] < metrics[(unsigned int)best_valid]))
                    best_valid = (int)b;
            }
            unsigned int teacher = best_valid >= 0 ? (unsigned int)best_valid : (unsigned int)best_any;

            // score every branch with the trace model, then rank by ITS
            // criterion (descending if direction_largest, i.e. basin-style
            // probability, ascending otherwise) -- replaces pm's
            // ascending-pm_min ranking.
            std::vector<double> scores(L);
            for (unsigned int b = 0; b < L; b++)
                scores[b] = scoreForBranch(model, b, checkpoint_trace[b]);
            std::vector<unsigned int> order(L);
            for (unsigned int b = 0; b < L; b++) order[b] = b;
            std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
                return largest ? scores[a] > scores[b] : scores[a] < scores[b];
            });

            // m_required = 1-indexed rank of `teacher` within the
            // trace-score order (same well-posedness argument as pm: a
            // global pickBest winner is the pickBest winner of any subset
            // containing it, so this is monotonic regardless of which
            // ranking criterion produced the prefix order).
            unsigned int m_required = 0;
            for (unsigned int r = 0; r < L; r++)
                if (order[r] == teacher) { m_required = r + 1; break; }

            out << s << "," << L << "," << m_required << "," << (correct[teacher] ? 1 : 0);
            for (unsigned int r = 0; r < L; r++) {
                double sc = scores[order[r]];
                double val = largest ? 1.0 / (1.0 + std::exp(-sc)) : sc;
                out << "," << val;
            }
            out << "\n";

            if ((s + 1) % 2000 == 0) {
                std::cerr << "[progress] samples=" << (s + 1) << "/" << args.samples << "\n";
                std::cerr.flush();
            }
        }
        out.close();
        std::cerr << "wrote " << args.out_path << " (samples=" << args.samples << ")\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
