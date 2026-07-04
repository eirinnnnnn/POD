#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "libParser.h"
#include "mclass_decoder_wrapper.h"

struct Args {
    std::string ini_path;
    std::string dataset_csv;
    std::string eps_csv = "1e-4,1e-3,1e-2,1e-1";
    unsigned int max_samples = 200;
    unsigned int trials = 5;
    unsigned int seed = 1;
};

struct Dataset {
    std::vector<std::string> header;
    std::map<std::string, unsigned int> col;
    std::vector<std::vector<std::string> > rows;
};

struct DecodeSummary {
    unsigned int label = 0;
    bool correct = false;
    std::vector<char> word;
    std::vector<char> basin;
};

struct Stats {
    unsigned long total = 0;
    unsigned long label_changed = 0;
    unsigned long word_changed = 0;
    unsigned long correctness_changed = 0;
    unsigned long pert_label_in_orig_basin = 0;
    unsigned long orig_label_in_pert_basin = 0;
    unsigned long basin_membership_changed = 0;
};

static std::vector<std::string> split(const std::string &line, char delim) {
    std::vector<std::string> out;
    std::string cur;
    std::stringstream ss(line);
    while(std::getline(ss, cur, delim))
        out.push_back(cur);
    return out;
}

static std::vector<double> parseDoubles(const std::string &s) {
    std::vector<double> out;
    for(const std::string &part : split(s, ','))
        if(!part.empty())
            out.push_back(std::stod(part));
    return out;
}

static Dataset readDataset(const std::string &path) {
    std::ifstream in(path.c_str());
    if(!in)
        throw std::runtime_error("cannot open dataset: " + path);
    Dataset data;
    std::string line;
    if(!std::getline(in, line))
        throw std::runtime_error("empty dataset: " + path);
    data.header = split(line, ',');
    for(unsigned int i=0; i<data.header.size(); i++)
        data.col[data.header[i]] = i;
    while(std::getline(in, line)) {
        if(!line.empty())
            data.rows.push_back(split(line, ','));
    }
    return data;
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

static Args parseArgs(int argc, char **argv) {
    Args args;
    for(int i=1; i<argc; i++) {
        std::string a = argv[i];
        auto need = [&](const std::string &name)->std::string {
            if(i + 1 >= argc)
                throw std::runtime_error("missing value after " + name);
            return argv[++i];
        };
        if(a == "-ini" || a == "--ini")
            args.ini_path = need(a);
        else if(a == "--dataset")
            args.dataset_csv = need(a);
        else if(a == "--eps")
            args.eps_csv = need(a);
        else if(a == "--max-samples")
            args.max_samples = (unsigned int)std::stoul(need(a));
        else if(a == "--trials")
            args.trials = (unsigned int)std::stoul(need(a));
        else if(a == "--seed")
            args.seed = (unsigned int)std::stoul(need(a));
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if(args.ini_path.empty())
        throw std::runtime_error("--ini is required");
    if(args.dataset_csv.empty())
        throw std::runtime_error("--dataset is required");
    return args;
}

static DecodeSummary decodeAll(MClassDecoder &decoder,
                               const std::vector<double> &received,
                               const std::vector<char> &codeword) {
    const unsigned int M = decoder.branchCount();
    const unsigned int n = decoder.codewordLength();
    std::string log;
    std::vector<double> metrics(M, std::numeric_limits<double>::max());
    std::vector<std::vector<char> > words(M, std::vector<char>(n, 0));
    std::vector<char> ok(M, 0);

    double best = std::numeric_limits<double>::max();
    unsigned int label = 0;
    for(unsigned int b=0; b<M; b++) {
        double metric = 0.0;
        ok[b] = decoder.decodeBranch(received, b, words[b], metric, log);
        metrics[b] = ok[b] ? metric : std::numeric_limits<double>::max();
        if(ok[b] && metric < best) {
            best = metric;
            label = b;
        }
    }

    DecodeSummary s;
    s.label = label;
    s.word = words[label];
    s.correct = (s.word == codeword);
    s.basin.assign(M, 0);
    for(unsigned int b=0; b<M; b++)
        if(ok[b] && std::abs(metrics[b] - best) <= 1e-12)
            s.basin[b] = 1;
    return s;
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);
        std::vector<double> eps_values = parseDoubles(args.eps_csv);
        Dataset data = readDataset(args.dataset_csv);

        std::map<std::string, std::map<std::string, std::string> > config = defaultConfig();
        parseConfig(args.ini_path, config);
        config["AdjustPolarDecoder"]["OnlyInit"] = "false";
        config["AdjustPolarDecoder"]["use_AED"] = "true";
        MClassDecoder decoder(config["AdjustPolarDecoder"]);
        const unsigned int n = decoder.codewordLength();
        const unsigned int rows = std::min<unsigned int>(args.max_samples, (unsigned int)data.rows.size());

        std::vector<Stats> stats(eps_values.size());
        std::mt19937 rng(args.seed);
        std::normal_distribution<double> normal(0.0, 1.0);

        for(unsigned int r=0; r<rows; r++) {
            const std::vector<std::string> &row = data.rows[r];
            std::vector<double> y(n, 0.0);
            std::vector<char> cw(n, 0);
            for(unsigned int i=0; i<n; i++) {
                y[i] = std::stod(row[data.col.at("y_" + std::to_string(i))]);
                cw[i] = (char)std::stoi(row[data.col.at("cw_" + std::to_string(i))]);
            }
            DecodeSummary orig = decodeAll(decoder, y, cw);

            for(unsigned int e=0; e<eps_values.size(); e++) {
                for(unsigned int t=0; t<args.trials; t++) {
                    std::vector<double> yp = y;
                    for(double &v : yp)
                        v += eps_values[e] * normal(rng);
                    DecodeSummary pert = decodeAll(decoder, yp, cw);

                    Stats &st = stats[e];
                    st.total++;
                    if(pert.label != orig.label)
                        st.label_changed++;
                    if(pert.word != orig.word)
                        st.word_changed++;
                    if(pert.correct != orig.correct)
                        st.correctness_changed++;
                    if(orig.basin[pert.label])
                        st.pert_label_in_orig_basin++;
                    if(pert.basin[orig.label])
                        st.orig_label_in_pert_basin++;
                    if(pert.basin != orig.basin)
                        st.basin_membership_changed++;
                }
            }
        }

        std::cout << "samples," << rows << "\n";
        std::cout << "trials," << args.trials << "\n";
        std::cout << "M," << decoder.branchCount() << "\n";
        std::cout << "eps,total,label_change,word_change,correctness_change,"
                     "pert_label_in_orig_basin,orig_label_in_pert_basin,basin_membership_change\n";
        for(unsigned int e=0; e<eps_values.size(); e++) {
            const Stats &st = stats[e];
            const double total = (double)st.total;
            std::cout << eps_values[e]
                      << "," << st.total
                      << "," << st.label_changed / total
                      << "," << st.word_changed / total
                      << "," << st.correctness_changed / total
                      << "," << st.pert_label_in_orig_basin / total
                      << "," << st.orig_label_in_pert_basin / total
                      << "," << st.basin_membership_changed / total
                      << "\n";
        }
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
