// Dataset generator for the "adaptive path selector, model_out variant":
// same SET-LEVEL target as adaptive_path_selector/pm (predict, per SAMPLE,
// the minimum pruning width k such that keeping the top-k branches still
// contains a branch full-PED/basin would have picked), but ranked and
// featurized by the trace_learnt_path_selector's OWN per-branch output
// instead of the raw checkpoint pm_min.
//
// trace_learnt_path_selector's model IO (see
// trace_learnt_path_selector/*/common/mclass_bler_sim.cpp): input per
// branch is that branch's path-metric AND channel-LLR-magnitude statistics
// (pm_min/pm_gap/pm_mean/pm_max, llr_abs_min/mean/max) at every checkpoint
// up to t, plus a branch-identity one-hot; output is a single scalar per
// branch (a BCEWithLogits logit for target=basin models -- sigmoid(logit)
// is then a literal per-branch correctness probability). Stacking that
// scalar over all M branches gives exactly the dim-M probability vector
// this generator uses as its OWN input feature, replacing pm's
// t{T}_pm_min_sorted_1..M columns with t{T}_trace_prob_sorted_1..M (sorted
// by the trace model's own ranking criterion, not by pm_min). T is taken
// from the trace model itself (its highest trained checkpoint), not a
// separate --checkpoint argument -- there is no independent ranking
// checkpoint here, only whatever checkpoints the loaded model was trained
// on.
//
// This is a NEW file only -- copies (does not modify) the Model/loadModel/
// scoreForBranch block from trace_learnt_path_selector's
// common/mclass_bler_sim.cpp, and reuses common/mclass_decoder_wrapper.h
// (MClassDecoder) unmodified. See adaptive_path_selector/pm/ebch_m7t10/
// mclass_adaptive_m_dataset.cpp for the pm-input sibling this mirrors.
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
// common/mclass_bler_sim.cpp (Model/loadModel/scoreForBranch) ----
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

// same defaults as adaptive_path_selector/pm/ebch_m7t10
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
        Model model = loadModel(args.trace_model_path);
        const unsigned int T = model.checkpoints_needed.back();

        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        const unsigned int n = decoder.codewordLength();
        const unsigned int k = decoder.messageLength();
        const unsigned int M = decoder.branchCount();
        if (T > n) throw std::runtime_error("trace model's max checkpoint exceeds codeword length");
        const std::vector<unsigned int> &checkpoints = model.checkpoints_needed;
        const bool largest = model.direction_largest;
        const char *col_base = largest ? "trace_prob_sorted" : "trace_score_sorted";

        long msg_seed = args.message_seed;
        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        std::ofstream out(args.out_path);
        if (!out) throw std::runtime_error("cannot open --out for writing: " + args.out_path);
        out << "sample,M,m_required,full_ped_correct,m_required_basin,any_basin_exists";
        for (unsigned int r = 1; r <= M; r++)
            out << ",t" << T << "_" << col_base << "_" << r;
        out << "\n";

        for (unsigned int s = 0; s < args.samples; s++) {
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

            // full_ped's own choice: global argmin metric, same convention
            // as adaptive_path_selector/pm and common/mclass_bler_sim.cpp
            // (tracked by VALUE -- exact metric ties share the same
            // correctness outcome, so which tied index is picked never
            // matters).
            unsigned int teacher = 0;
            for (unsigned int b = 1; b < M; b++)
                if (metrics[b] < metrics[teacher]) teacher = b;
            const double global_min_metric = metrics[teacher];

            // score every branch with the trace model, then rank by ITS
            // criterion (descending if direction_largest, i.e. basin-style
            // probability, ascending otherwise) -- this replaces pm's
            // ascending-pm_min ranking.
            std::vector<double> scores(M);
            for (unsigned int b = 0; b < M; b++)
                scores[b] = scoreForBranch(model, b, traces[b]);
            std::vector<unsigned int> order(M);
            for (unsigned int b = 0; b < M; b++) order[b] = b;
            std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
                return largest ? scores[a] > scores[b] : scores[a] < scores[b];
            });

            // m_required / m_required_basin: same definitions as pm, just
            // evaluated against the trace-score order instead of pm_min order.
            unsigned int m_required = 0;
            double running_min = std::numeric_limits<double>::max();
            for (unsigned int r = 0; r < M; r++) {
                if (metrics[order[r]] < running_min) running_min = metrics[order[r]];
                if (running_min <= global_min_metric) { m_required = r + 1; break; }
            }

            unsigned int m_required_basin = 0;
            bool any_basin_exists = false;
            for (unsigned int r = 0; r < M; r++) {
                if (correct[order[r]]) { m_required_basin = r + 1; break; }
            }
            for (unsigned int b = 0; b < M; b++)
                if (correct[b]) { any_basin_exists = true; break; }

            out << s << "," << M << "," << m_required << "," << (correct[teacher] ? 1 : 0)
                << "," << m_required_basin << "," << (any_basin_exists ? 1 : 0);
            for (unsigned int r = 0; r < M; r++) {
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
