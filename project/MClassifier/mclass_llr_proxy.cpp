#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "libParser.h"
#include "mclass_decoder_wrapper.h"

struct Args {
    std::string ini_path;
    std::string dataset_csv;
    std::string out_csv;
    std::string selection = "metric_argmin";
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
    config["AdjustPolarDecoder"]["use_AED"] = "true";
    config["AdjustPolarDecoder"]["automorphism_src"] = "";
    config["AdjustPolarDecoder"]["aed_L"] = "0";

    return config;
}

static void usage() {
    std::cout
        << "Usage:\n"
        << "  mclass_llr_proxy --ini config.ini --dataset dataset.csv --out proxy.csv [options]\n\n"
        << "Options:\n"
        << "  --selection FIELD       metric_argmin, valid_metric_argmin, or label. Default metric_argmin.\n";
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
        } else if(a == "--dataset") {
            args.dataset_csv = need(a);
        } else if(a == "--out") {
            args.out_csv = need(a);
        } else if(a == "--selection") {
            args.selection = need(a);
        } else if(a == "-h" || a == "--help") {
            usage();
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + a);
        }
    }

    if(args.ini_path.empty())
        throw std::runtime_error("--ini is required");
    if(args.dataset_csv.empty())
        throw std::runtime_error("--dataset is required");
    if(args.out_csv.empty())
        throw std::runtime_error("--out is required");
    return args;
}

static std::vector<std::string> splitCsvLine(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while(std::getline(ss, item, ','))
        out.push_back(item);
    return out;
}

static int findColumn(const std::vector<std::string> &header, const std::string &name) {
    for(unsigned int i=0; i<header.size(); i++)
        if(header[i] == name)
            return (int)i;
    return -1;
}

static double scoreBranch(const std::vector<double> &y,
                          const std::vector<unsigned int> &order,
                          const std::vector<unsigned int> &info_pos,
                          unsigned int prefix) {
    std::vector<unsigned int> src_for_pos(order.size(), 0);
    for(unsigned int src=0; src<order.size(); src++)
        src_for_pos[order[src]] = src;

    unsigned int count = std::min<unsigned int>(prefix, (unsigned int)info_pos.size());
    double sum = 0.0;
    for(unsigned int i=0; i<count; i++)
        sum += std::abs(y[src_for_pos[info_pos[i]]]);
    return count ? sum / count : 0.0;
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);

        std::map<std::string, std::map<std::string, std::string> > config = defaultConfig();
        parseConfig(args.ini_path, config);
        config["AdjustPolarDecoder"]["OnlyInit"] = "false";
        config["AdjustPolarDecoder"]["use_AED"] = "true";

        MClassDecoder decoder(config["AdjustPolarDecoder"]);
        const unsigned int n = decoder.codewordLength();
        const unsigned int M = decoder.branchCount();
        std::vector<unsigned int> info_pos = decoder.informationPositions();
        std::vector<unsigned int> prefixes;
        unsigned int candidates[] = {8, 16, 32, 64};
        for(unsigned int i=0; i<4; i++)
            if(candidates[i] <= info_pos.size())
                prefixes.push_back(candidates[i]);
        if(prefixes.empty())
            prefixes.push_back((unsigned int)info_pos.size());

        std::ifstream in(args.dataset_csv.c_str());
        if(!in)
            throw std::runtime_error("cannot open dataset: " + args.dataset_csv);
        std::ofstream out(args.out_csv.c_str());
        if(!out)
            throw std::runtime_error("cannot open output: " + args.out_csv);
        out << std::setprecision(17);

        std::string line;
        if(!std::getline(in, line))
            throw std::runtime_error("empty dataset");
        std::vector<std::string> header = splitCsvLine(line);

        int sample_col = findColumn(header, "sample");
        int select_col = findColumn(header, args.selection);
        int teacher_col = findColumn(header, "teacher_correct");
        int y0_col = findColumn(header, "y_0");
        if(sample_col < 0 || select_col < 0 || y0_col < 0)
            throw std::runtime_error("dataset missing sample/selection/y_0 columns");

        out << "sample,selected_branch";
        if(teacher_col >= 0)
            out << ",teacher_correct";
        for(unsigned int p=0; p<prefixes.size(); p++) {
            out << ",sel_score_info" << prefixes[p]
                << ",best_score_info" << prefixes[p]
                << ",avg_score_info" << prefixes[p]
                << ",sel_rank_info" << prefixes[p]
                << ",best_branch_info" << prefixes[p];
        }
        out << "\n";

        unsigned int rows = 0;
        while(std::getline(in, line)) {
            if(line.empty())
                continue;
            std::vector<std::string> row = splitCsvLine(line);
            if(row.size() != header.size())
                throw std::runtime_error("malformed csv row");

            std::vector<double> y(n, 0.0);
            for(unsigned int i=0; i<n; i++)
                y[i] = std::stod(row[y0_col + (int)i]);

            unsigned int selected = (unsigned int)std::stoul(row[select_col]);
            out << row[sample_col] << "," << selected;
            if(teacher_col >= 0)
                out << "," << row[teacher_col];

            for(unsigned int p=0; p<prefixes.size(); p++) {
                std::vector<double> scores(M, 0.0);
                double sum = 0.0;
                unsigned int best_branch = 0;
                for(unsigned int b=0; b<M; b++) {
                    scores[b] = scoreBranch(y, decoder.branchOrder(b), info_pos, prefixes[p]);
                    sum += scores[b];
                    if(scores[b] > scores[best_branch])
                        best_branch = b;
                }
                unsigned int rank = 1;
                for(unsigned int b=0; b<M; b++) {
                    if(scores[b] > scores[selected])
                        rank++;
                }
                out << "," << scores[selected]
                    << "," << scores[best_branch]
                    << "," << (sum / M)
                    << "," << rank
                    << "," << best_branch;
            }
            out << "\n";
            rows++;
        }

        std::cout << "wrote " << args.out_csv << "\n";
        std::cout << "rows=" << rows << " n=" << n << " k=" << info_pos.size()
                  << " M=" << M << " selection=" << args.selection << "\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
