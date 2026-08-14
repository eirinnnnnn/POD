// Trace-feature dataset generator for the Golay24 trace-learned branch
// ranker. For each sample, decodes all M AED branches with a LIST of
// checkpoint snapshots (cumulative -- e.g. --checkpoints 4,8,16 captures
// all three, matching the eBCH "upto_tX" convention where a t=T model sees
// every checkpoint <= T, not just T itself) and emits one dense row per
// (sample, branch) with each checkpoint's features plus training targets,
// matching the column convention expected by
// architectures/trace/train_mclass_trace_ranker.py:
//   sample, branch,
//   t{c}_active, t{c}_pm_min, t{c}_pm_gap, t{c}_pm_mean, t{c}_pm_max,
//   t{c}_llr_abs_min, t{c}_llr_abs_mean, t{c}_llr_abs_max   (repeated per checkpoint c)
//   metric_gap, basin, correct, final_metric
// basin = correct (1.0 if this branch's own decode equals the transmitted
// codeword) -- the natural per-branch label for "would keeping this branch
// let us recover the right answer".
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
    std::vector<unsigned int> checkpoints;
    long channel_seed = -2;
    long message_seed = -1;
    std::string out_path;
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
        else if (s == "--checkpoint") a.checkpoints = {(unsigned int)std::stoul(need(s))};
        else if (s == "--checkpoints") a.checkpoints = parseUIntList(need(s));
        else if (s == "--channel-seed") a.channel_seed = std::stol(need(s));
        else if (s == "--message-seed") a.message_seed = std::stol(need(s));
        else if (s == "--out") a.out_path = need(s);
        else throw std::runtime_error("unknown argument: " + s);
    }
    if (a.ini_path.empty()) throw std::runtime_error("-ini is required");
    if (!a.snr_set) throw std::runtime_error("--snr is required");
    if (a.samples == 0) throw std::runtime_error("--samples is required");
    if (a.checkpoints.empty()) throw std::runtime_error("--checkpoints (or --checkpoint) is required");
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
        for (unsigned int c : args.checkpoints)
            if (c > n) throw std::runtime_error("--checkpoints contains a value exceeding codeword length");
        const std::vector<unsigned int> &checkpoints = args.checkpoints;
        const unsigned int nc = (unsigned int)checkpoints.size();

        long msg_seed = args.message_seed;
        std::vector<char> message(k, 0), codeword(n, 0);
        std::vector<double> received(n, 0.0);

        std::ofstream out(args.out_path);
        if (!out) throw std::runtime_error("cannot open --out for writing: " + args.out_path);
        out << "sample,branch";
        for (unsigned int c : checkpoints)
            out << ",t" << c << "_active,t" << c << "_pm_min,t" << c << "_pm_gap,"
                << "t" << c << "_pm_mean,t" << c << "_pm_max,"
                << "t" << c << "_llr_abs_min,t" << c << "_llr_abs_mean,t" << c << "_llr_abs_max";
        out << ",metric_gap,basin,correct,final_metric\n";

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
            const unsigned int M = (unsigned int)candidates.size();

            double min_metric = 1e100;
            for (unsigned int b = 0; b < M; b++)
                if (metrics[b] < min_metric) min_metric = metrics[b];

            for (unsigned int b = 0; b < M; b++) {
                bool correct = (candidates[b] == codeword);
                double metric_gap = metrics[b] - min_metric;
                out << s << "," << b;
                for (unsigned int ci = 0; ci < nc; ci++) {
                    const PSPEDTracePoint &tp = checkpoint_trace[b][ci];
                    out << "," << tp.active << "," << tp.pm_min << "," << tp.pm_gap << ","
                        << tp.pm_mean << "," << tp.pm_max << ","
                        << tp.llr_abs_min << "," << tp.llr_abs_mean << "," << tp.llr_abs_max;
                }
                out << "," << metric_gap << "," << (correct ? 1 : 0) << "," << (correct ? 1 : 0) << ","
                    << metrics[b] << "\n";
            }

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
