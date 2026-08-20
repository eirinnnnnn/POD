// Dataset generator for the "adaptive path selector": a SET-LEVEL model,
// stacked on top of the static-pm ranker, that predicts -- per sample, not
// per branch -- the minimum pruning width m such that keeping the top-m
// branches (sorted ascending by their CHECKPOINT-t pm_min, exactly the
// static-pm-rank criterion in golay_mclass_bler_sim.cpp) still contains the
// branch that full-PED itself would have picked.
//
// "Teacher" = pickBest(all branches): prefer a parity-valid candidate with
// the lowest FINAL metric, else fall back to the lowest final metric
// overall. This is byte-for-byte the same policy golay_mclass_bler_sim.cpp
// uses for its full_ped reference curve and for pruned-subset selection --
// reusing any other definition (e.g. "the branch whose own decode is
// correct") would silently diverge from what static-pm-rank/trace-learned
// are actually being compared against.
//
// One row per SAMPLE (not per branch):
//   sample, m_required, L, teacher_is_correct, full_ped_correct,
//   t{T}_pm_min_sorted_1 .. t{T}_pm_min_sorted_L
// The sorted pm_min vector at checkpoint T is the adaptive selector's
// input: it deliberately drops which physical branch is which, since
// m_required is defined purely by the SHAPE of that sorted score profile
// (how sharply the front-runner separates from the pack) -- it is a
// sufficient statistic by construction.
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
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
    unsigned int checkpoint = 0;   // single ranking checkpoint t (decode_idx progress count)
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
        else if (s == "--checkpoint") a.checkpoint = (unsigned int)std::stoul(need(s));
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--out") a.out_path = need(s);
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (!a.snr_set) throw std::runtime_error("--snr is required");
    if (a.samples == 0) throw std::runtime_error("--samples is required");
    if (a.checkpoint == 0) throw std::runtime_error("--checkpoint is required (t, decode_idx progress count)");
    if (a.out_path.empty()) throw std::runtime_error("--out is required");
    return a;
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
        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        const unsigned int n = decoder.getCodewordLength();
        const unsigned int k = decoder.getMessageLength();
        if (args.checkpoint > n) throw std::runtime_error("--checkpoint exceeds codeword length");
        std::vector<unsigned int> checkpoints;
        { static const unsigned int ladder[] = {4,8,16,20};
          for (unsigned int c : ladder) if (c <= args.checkpoint) checkpoints.push_back(c);
          if (checkpoints.empty() || checkpoints.back() != args.checkpoint) checkpoints.push_back(args.checkpoint); }
        const unsigned int rank_cp_pos = (unsigned int)(std::lower_bound(checkpoints.begin(), checkpoints.end(), args.checkpoint) - checkpoints.begin());

        long msg_seed = args.message_seed;
        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        const unsigned int L_fixed = decoder.branchCount();

        std::ofstream out(args.out_path);
        if (!out) throw std::runtime_error("cannot open --out for writing: " + args.out_path);
        out << "sample,L,m_required,full_ped_correct";
        for (unsigned int r = 1; r <= L_fixed; r++)
            out << ",t" << args.checkpoint << "_pm_min_sorted_" << r;
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

            // teacher = pickBest(all L branches), byte-for-byte the same
            // policy golay_mclass_bler_sim.cpp uses for its full_ped
            // reference: prefer a parity-valid candidate with the lowest
            // FINAL metric, else lowest final metric overall.
            int best_valid = -1, best_any = -1;
            for (unsigned int b = 0; b < L; b++) {
                if (best_any < 0 || metrics[b] < metrics[(unsigned int)best_any]) best_any = (int)b;
                if (valid[b] && (best_valid < 0 || metrics[b] < metrics[(unsigned int)best_valid]))
                    best_valid = (int)b;
            }
            unsigned int teacher = best_valid >= 0 ? (unsigned int)best_valid : (unsigned int)best_any;

            // sort branch indices ascending by checkpoint pm_min (the
            // static-pm-rank criterion)
            std::vector<unsigned int> order(L);
            for (unsigned int b = 0; b < L; b++) order[b] = b;
            std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
                return checkpoint_trace[a][rank_cp_pos].pm_min < checkpoint_trace[b][rank_cp_pos].pm_min;
            });

            // m_required = 1-indexed rank of `teacher` within that sorted
            // order, i.e. the smallest m for which keeping the top-m
            // branches (by checkpoint pm_min) is GUARANTEED to reproduce
            // full_ped's own decision exactly. This is well-posed and
            // monotonic: teacher is, by construction, the branch with the
            // globally lowest metric among its own preference class (valid
            // if any exist, else any) across ALL L branches, so once it
            // enters a growing prefix it remains that prefix's pickBest()
            // result for every larger prefix too (a global min is the
            // local min of any subset containing it). Note this is NOT the
            // same quantity as static_pm_rank(t,m)'s raw per-m correctness
            // in the BLER sweeps: raw correctness at a FIXED m can flicker
            // non-monotonically as m grows (a wrong-but-parity-valid
            // competitor with an even lower metric can knock out a
            // correct answer once more branches are admitted), so it is
            // not a well-posed single-threshold target. m_required instead
            // answers "how many branches, sorted by checkpoint pm_min, are
            // enough to exactly match what full-PED itself would decide"
            // -- which is the quantity that actually lets an adaptive-m
            // scheme reproduce full_ped's BLER exactly (not just
            // approximate it) whenever the prediction is accurate.
            unsigned int m_required = 0;
            for (unsigned int r = 0; r < L; r++)
                if (order[r] == teacher) { m_required = r + 1; break; }

            out << s << "," << L << "," << m_required << "," << (correct[teacher] ? 1 : 0);
            for (unsigned int r = 0; r < L; r++)
                out << "," << checkpoint_trace[order[r]][rank_cp_pos].pm_min;
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
