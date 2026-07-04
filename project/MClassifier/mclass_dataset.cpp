#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
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
    std::string out_csv;
    std::string label_mode = "oracle_correct_then_metric";
    unsigned int samples = 0;
    std::string snr;
    std::string channel_seed;
    std::string message_seed;
};

struct BranchResult {
    bool ok = false;
    bool valid = false;
    bool correct = false;
    double metric = 0.0;
    std::vector<char> word;
};

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
    config["AdjustPolarDecoder"]["Gmatrix_path"] = "../../../data/worse_case_8bit.matrix";
    config["AdjustPolarDecoder"]["permutation_src"] = "";
    config["AdjustPolarDecoder"]["permutation_path"] = "";
    config["AdjustPolarDecoder"]["permutation_list"] = "";
    config["AdjustPolarDecoder"]["permutation_random_seed"] = "-1";
    config["AdjustPolarDecoder"]["target_raw_BER"] = "0.01";
    config["AdjustPolarDecoder"]["operationArray"] = "111011111111";
    config["AdjustPolarDecoder"]["list_size"] = "4";
    config["AdjustPolarDecoder"]["bha_value_setting"] = "????????";
    config["AdjustPolarDecoder"]["OnlyInit"] = "false";
    config["AdjustPolarDecoder"]["use_AED"] = "false";
    config["AdjustPolarDecoder"]["automorphism_src"] = "";
    config["AdjustPolarDecoder"]["aed_L"] = "0";

    return config;
}

static void usage() {
    std::cout
        << "Usage:\n"
        << "  mclass_dataset -ini config.ini --samples N --out dataset.csv [options]\n\n"
        << "Options:\n"
        << "  --snr VALUE                  Override AWGN start/end to one SNR.\n"
        << "  --channel-seed VALUE         Override AWGN seed_string.\n"
        << "  --message-seed VALUE         Override Monte_Carlo seed_string.\n"
        << "  --label-mode MODE            oracle_correct_then_metric, metric, valid_metric.\n";
}

static Args parseArgs(int argc, char **argv) {
    Args args;
    for(int i=1; i<argc; i++) {
        std::string a = argv[i];
        auto need = [&](const std::string &name)->std::string {
            if(i + 1 >= argc)
                throw std::runtime_error("missing value after " + name);
            return argv[++i];
        };

        if(a == "-ini" || a == "--ini") {
            args.ini_path = need(a);
        } else if(a == "--samples") {
            args.samples = (unsigned int)std::stoul(need(a));
        } else if(a == "--out") {
            args.out_csv = need(a);
        } else if(a == "--snr") {
            args.snr = need(a);
        } else if(a == "--channel-seed") {
            args.channel_seed = need(a);
        } else if(a == "--message-seed") {
            args.message_seed = need(a);
        } else if(a == "--label-mode") {
            args.label_mode = need(a);
        } else if(a == "-h" || a == "--help") {
            usage();
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + a);
        }
    }

    if(args.ini_path.empty())
        throw std::runtime_error("-ini is required");
    if(args.out_csv.empty())
        throw std::runtime_error("--out is required");
    if(args.samples == 0)
        throw std::runtime_error("--samples must be positive");
    return args;
}

static std::string metadataPath(const std::string &csv_path) {
    return csv_path + ".meta.json";
}

static unsigned int argminMetric(const std::vector<BranchResult> &branches,
                                 bool require_valid,
                                 bool require_correct,
                                 bool *found) {
    double best = std::numeric_limits<double>::max();
    unsigned int best_idx = 0;
    bool have = false;

    for(unsigned int i=0; i<branches.size(); i++) {
        const BranchResult &b = branches[i];
        if(!b.ok)
            continue;
        if(require_valid && !b.valid)
            continue;
        if(require_correct && !b.correct)
            continue;
        if(!have || b.metric < best) {
            have = true;
            best = b.metric;
            best_idx = i;
        }
    }

    if(found)
        *found = have;
    return best_idx;
}

static unsigned int chooseLabel(const std::vector<BranchResult> &branches,
                                const std::string &mode) {
    bool found = false;
    if(mode == "metric") {
        return argminMetric(branches, false, false, &found);
    }
    if(mode == "valid_metric") {
        unsigned int idx = argminMetric(branches, true, false, &found);
        if(found)
            return idx;
        return argminMetric(branches, false, false, &found);
    }
    if(mode == "oracle_correct_then_metric") {
        unsigned int idx = argminMetric(branches, false, true, &found);
        if(found)
            return idx;
        return argminMetric(branches, false, false, &found);
    }
    throw std::runtime_error("unknown label mode: " + mode);
}

static void writeMetadata(const Args &args,
                          const std::map<std::string, std::map<std::string, std::string> > &config,
                          MClassDecoder &decoder) {
    std::ofstream meta(metadataPath(args.out_csv).c_str());
    if(!meta)
        throw std::runtime_error("cannot open metadata output");

    meta << "{\n";
    meta << "  \"format\": \"mclass_csv_v1\",\n";
    meta << "  \"csv_path\": \"" << args.out_csv << "\",\n";
    meta << "  \"ini_path\": \"" << args.ini_path << "\",\n";
    meta << "  \"samples\": " << args.samples << ",\n";
    meta << "  \"label_mode\": \"" << args.label_mode << "\",\n";
    meta << "  \"n\": " << decoder.codewordLength() << ",\n";
    meta << "  \"k\": " << decoder.messageLength() << ",\n";
    meta << "  \"M\": " << decoder.branchCount() << ",\n";
    meta << "  \"list_size\": " << decoder.sclListSize() << ",\n";
    meta << "  \"awgn_start\": \"" << config.at("AWGN").at("start") << "\",\n";
    meta << "  \"awgn_step_type\": \"" << config.at("AWGN").at("step_type") << "\",\n";
    meta << "  \"awgn_seed\": \"" << config.at("AWGN").at("seed_string") << "\",\n";
    meta << "  \"message_seed\": \"" << config.at("Monte_Carlo").at("seed_string") << "\"\n";
    meta << "}\n";
}

static void writeHeader(std::ofstream &out, unsigned int n, unsigned int k, unsigned int M) {
    out << "sample,label,metric_argmin,valid_metric_argmin,teacher_correct";
    for(unsigned int i=0; i<k; i++)
        out << ",msg_" << i;
    for(unsigned int i=0; i<n; i++)
        out << ",cw_" << i;
    for(unsigned int i=0; i<n; i++)
        out << ",y_" << i;
    for(unsigned int i=0; i<M; i++)
        out << ",metric_" << i;
    for(unsigned int i=0; i<M; i++)
        out << ",valid_" << i;
    for(unsigned int i=0; i<M; i++)
        out << ",correct_" << i;
    out << "\n";
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);

        std::map<std::string, std::map<std::string, std::string> > config = defaultConfig();
        parseConfig(args.ini_path, config);
        config["AdjustPolarDecoder"]["OnlyInit"] = "false";
        config["AdjustPolarDecoder"]["use_AED"] = "true";
        if(!args.snr.empty()) {
            config["AWGN"]["step_type"] = "SNR";
            config["AWGN"]["start"] = args.snr;
            config["AWGN"]["end"] = args.snr;
            config["AWGN"]["step"] = "1.0";
        }
        if(!args.channel_seed.empty())
            config["AWGN"]["seed_string"] = args.channel_seed;
        if(!args.message_seed.empty())
            config["Monte_Carlo"]["seed_string"] = args.message_seed;

        MClassDecoder decoder(config["AdjustPolarDecoder"]);
        if(decoder.branchCount() == 0)
            throw std::runtime_error("decoder has zero PED branches");

        AWGN channel(config["AWGN"]);
        channel.setCodeRate(decoder.getCodeRate());

        const unsigned int k = decoder.messageLength();
        const unsigned int n = decoder.codewordLength();
        const unsigned int M = decoder.branchCount();

        std::ofstream out(args.out_csv.c_str());
        if(!out)
            throw std::runtime_error("cannot open output csv: " + args.out_csv);
        out << std::setprecision(17);
        writeHeader(out, n, k, M);

        long msg_seed = std::stol(config["Monte_Carlo"]["seed_string"]);
        std::vector<char> message(k, 0);
        std::vector<char> codeword(n, 0);
        std::vector<double> received(n, 0.0);
        std::string log;

        for(unsigned int sample=0; sample<args.samples; sample++) {
            for(unsigned int i=0; i<k; i++)
                message[i] = (ran0(&msg_seed) > 0.5 ? 1 : 0);

            decoder.doEncode(message, codeword);
            channel.addNoise(codeword, received);

            std::vector<BranchResult> branches(M);
            for(unsigned int b=0; b<M; b++) {
                branches[b].word.assign(n, 0);
                branches[b].ok = decoder.decodeBranch(received, b, branches[b].word, branches[b].metric, log);
                branches[b].valid = branches[b].ok && decoder.parityValid(branches[b].word);
                branches[b].correct = branches[b].ok && (branches[b].word == codeword);
            }

            bool found = false;
            unsigned int metric_argmin = argminMetric(branches, false, false, &found);
            unsigned int valid_metric_argmin = argminMetric(branches, true, false, &found);
            if(!found)
                valid_metric_argmin = metric_argmin;
            unsigned int label = chooseLabel(branches, args.label_mode);
            bool teacher_correct = branches[label].ok && branches[label].correct;

            out << sample
                << "," << label
                << "," << metric_argmin
                << "," << valid_metric_argmin
                << "," << (teacher_correct ? 1 : 0);
            for(unsigned int i=0; i<k; i++)
                out << "," << (int)message[i];
            for(unsigned int i=0; i<n; i++)
                out << "," << (int)codeword[i];
            for(unsigned int i=0; i<n; i++)
                out << "," << received[i];
            for(unsigned int i=0; i<M; i++)
                out << "," << (branches[i].ok ? branches[i].metric : 1e300);
            for(unsigned int i=0; i<M; i++)
                out << "," << (branches[i].valid ? 1 : 0);
            for(unsigned int i=0; i<M; i++)
                out << "," << (branches[i].correct ? 1 : 0);
            out << "\n";
        }

        writeMetadata(args, config, decoder);
        std::cout << "wrote " << args.out_csv << "\n";
        std::cout << "wrote " << metadataPath(args.out_csv) << "\n";
        std::cout << "samples=" << args.samples << " n=" << n << " k=" << k << " M=" << M << "\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
