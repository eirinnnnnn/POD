// Dataset generator for the "adaptive path selector": a SET-LEVEL model,
// stacked on top of the static-pm ranker (architectures/static_pm_rank),
// that predicts -- per SAMPLE, not per branch -- the minimum pruning width
// k such that keeping the top-k branches (sorted ascending by their
// CHECKPOINT-t pm_min, exactly the static-pm-rank criterion used in
// common/mclass_bler_sim.cpp) still contains the branch that full-PED
// itself would have picked.
//
// This is a NEW file only -- it reuses common/mclass_decoder_wrapper.h
// (MClassDecoder, already a non-invasive subclass of the original
// AdjustPolarDecoderRelation) unmodified, and does not touch any original
// decoder source.
//
// "Teacher" = full_ped's own choice = argmin(metrics) over all M branches
// -- byte-for-byte the same policy common/mclass_bler_sim.cpp uses for its
// full_ped reference AND for its static-pm-rank prefix scan (see the `best`
// ratchet in its per-k loop: no parity-valid preference is used anywhere
// in this codebase's eBCH tools, unlike Golay24's PS_PED-based sim, which
// does prefer parity-valid candidates -- this file matches eBCH's own
// convention, not Golay24's).
//
// One row per SAMPLE (not per branch):
//   sample, M, m_required, full_ped_correct, m_required_basin, any_basin_exists
//   t{T}_pm_min_sorted_1 .. t{T}_pm_min_sorted_M
//
// m_required = 1-indexed rank of the teacher branch within the checkpoint
// pm_min ascending sort. This is well-posed and MONOTONIC: teacher is, by
// construction, the branch with the globally lowest metric among ALL M
// branches, so once it enters a growing prefix (sorted by checkpoint
// pm_min) it remains that prefix's argmin for every larger prefix too (a
// global min is the local min of any subset containing it).
//
// m_required_basin = 1-indexed rank of the FIRST branch (by checkpoint
// pm_min order) whose OWN decode is correct -- not necessarily the specific
// metric-argmin branch full_ped would pick. Multiple AED branches can
// independently converge to the same correct codeword, so this is the
// more directly useful, generally SMALLER target: we don't need the
// pruned survivors to reproduce full_ped's exact metric-optimal choice,
// just to contain SOME correct candidate. any_basin_exists is 0 only when
// no branch among all M ever decodes correctly (a harder failure mode than
// full_ped merely picking the wrong one of several available correct
// branches). See NOTES.md for the full derivation and why m_required alone
// (matching a specific teacher branch) turned out to be an unnecessarily
// strict target.
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
#include "mclass_decoder_wrapper.h"

struct Args {
    std::string ini_path;
    double snr = 0.0;
    bool snr_set = false;
    unsigned int samples = 0;      // hard cap, always respected
    unsigned int target_errors = 0;  // 0 = disabled; else stop early once this many
                                      // full_ped_correct==0 events are collected
    unsigned int checkpoint = 0;   // single ranking checkpoint t (decode_idx progress count, 1..codeword_length)
    long channel_seed = -2;
    long message_seed = -1;
    std::string out_path;
    bool dump_frozen = false;   // diagnostic: print diverge_flag (info vs frozen per decode_idx) + relation span, and exit
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
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--out") a.out_path = need(s);
        else if (s == "--dump-frozen") a.dump_frozen = true;
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (!a.snr_set) throw std::runtime_error("--snr is required");
    if (a.samples == 0 && !a.dump_frozen) throw std::runtime_error("--samples is required");
    if (a.checkpoint == 0 && !a.dump_frozen) throw std::runtime_error("--checkpoint is required (t, decode_idx progress count)");
    if (a.out_path.empty() && !a.dump_frozen) throw std::runtime_error("--out is required");
    return a;
}

// same defaults as common/mclass_bler_sim.cpp
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

        if (args.dump_frozen) {
            std::vector<unsigned int> spans = decoder.relationSizes();
            std::vector<unsigned int> info_positions = decoder.informationPositions();
            std::vector<char> is_info(spans.size(), 0);
            for (unsigned int p : info_positions) is_info[p] = 1;
            for (unsigned int i = 0; i < spans.size(); i++) {
                std::cout << "decode_idx=" << i << " " << (is_info[i] ? "info" : "frozen")
                          << " relation_span=" << spans[i] << " relation=[";
                const std::vector<unsigned int> &rel = decoder.relationList(i);
                for (unsigned int j = 0; j < rel.size(); j++)
                    std::cout << (j ? "," : "") << rel[j];
                std::cout << "]\n";
            }
            return 0;
        }

        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        const unsigned int n = decoder.codewordLength();
        const unsigned int k = decoder.messageLength();
        const unsigned int M = decoder.branchCount();
        if (args.checkpoint > n) throw std::runtime_error("--checkpoint exceeds codeword length");
        std::vector<unsigned int> checkpoints;
        { static const unsigned int ladder[] = {8,16,32,64,128};
          for (unsigned int c : ladder) if (c <= args.checkpoint) checkpoints.push_back(c);
          if (checkpoints.empty() || checkpoints.back() != args.checkpoint) checkpoints.push_back(args.checkpoint); }
        const unsigned int rank_cp_pos = (unsigned int)(std::lower_bound(checkpoints.begin(), checkpoints.end(), args.checkpoint) - checkpoints.begin());

        long msg_seed = args.message_seed;
        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        std::ofstream out(args.out_path);
        if (!out) throw std::runtime_error("cannot open --out for writing: " + args.out_path);
        out << "sample,M,m_required,full_ped_correct,m_required_basin,any_basin_exists";
        for (unsigned int r = 1; r <= M; r++)
            out << ",t" << args.checkpoint << "_pm_min_sorted_" << r;
        out << "\n";

        unsigned int error_count = 0;
        unsigned int s = 0;
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

            // global_min_metric = full_ped's own metric = min(metrics) over
            // all M, no parity-valid preference (matches mclass_bler_sim.cpp).
            // Tracked by VALUE, not by a single "teacher" index: exact
            // metric ties are common (independently-converging AED branches
            // that reach the same correct codeword report bit-identical
            // metrics), and every branch tied for the global minimum shares
            // the same `correct` outcome, so which specific tied index is
            // picked never changes correctness -- only the achieved metric
            // VALUE does.
            unsigned int teacher = 0;
            for (unsigned int b = 1; b < M; b++)
                if (metrics[b] < metrics[teacher]) teacher = b;
            const double global_min_metric = metrics[teacher];

            // sort branch indices ascending by checkpoint pm_min (the
            // static-pm-rank criterion; single requested checkpoint, so
            // index 0 is t=args.checkpoint)
            std::vector<unsigned int> order(M);
            for (unsigned int b = 0; b < M; b++) order[b] = b;
            std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
                return traces[a][rank_cp_pos].pm_min < traces[b][rank_cp_pos].pm_min;
            });

            // m_required = smallest r such that the top-r branches (by
            // checkpoint pm_min) include AT LEAST ONE branch achieving the
            // global minimum metric -- equivalently, the smallest r at
            // which static-pm-rank(t,k>=r)'s running best-by-metric ratchet
            // locks onto full_ped's own metric value and, by the tie
            // argument above, its correctness outcome, for every k>=r.
            unsigned int m_required = 0;
            double running_min = std::numeric_limits<double>::max();
            for (unsigned int r = 0; r < M; r++) {
                if (metrics[order[r]] < running_min) running_min = metrics[order[r]];
                if (running_min <= global_min_metric) { m_required = r + 1; break; }
            }

            // m_required_basin = smallest r such that the top-r branches
            // (by checkpoint pm_min) contain ANY branch whose OWN decode is
            // correct (basin=1), not specifically full_ped's metric-argmin
            // choice. Multiple AED branches can independently converge to
            // the same correct codeword, so this can be much smaller than
            // m_required: we don't need to wait for the metric-optimal
            // branch specifically, just for any correct one to survive the
            // pruning. any_basin_exists (0 if no branch among all M ever
            // decodes correctly) is the natural companion diagnostic to
            // full_ped_correct=correct[teacher] -- they differ whenever
            // full_ped's own metric-argmin pick is wrong but some other,
            // non-metric-optimal branch was actually correct.
            unsigned int m_required_basin = 0;
            bool any_basin_exists = false;
            for (unsigned int r = 0; r < M; r++) {
                if (correct[order[r]]) { m_required_basin = r + 1; break; }
            }
            for (unsigned int b = 0; b < M; b++)
                if (correct[b]) { any_basin_exists = true; break; }

            out << s << "," << M << "," << m_required << "," << (correct[teacher] ? 1 : 0)
                << "," << m_required_basin << "," << (any_basin_exists ? 1 : 0);
            for (unsigned int r = 0; r < M; r++)
                out << "," << traces[order[r]][rank_cp_pos].pm_min;
            out << "\n";

            if (!correct[teacher]) error_count++;

            if ((s + 1) % 2000 == 0) {
                std::cerr << "[progress] samples=" << (s + 1) << "/" << args.samples
                           << " errors=" << error_count << "\n";
                std::cerr.flush();
            }
            if (args.target_errors > 0 && error_count >= args.target_errors) {
                s++;  // count this sample in the final total
                break;
            }
        }
        out.close();
        std::cerr << "wrote " << args.out_path << " (samples=" << s << ", errors=" << error_count << ")\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
