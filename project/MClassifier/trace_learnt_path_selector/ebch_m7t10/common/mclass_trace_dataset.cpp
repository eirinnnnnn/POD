#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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
    std::string checkpoints_csv = "8,16,32,64,128";
    unsigned int max_samples = 0;
};

// Dataset rows are read one line at a time (see StreamingDataset below)
// instead of being materialized into a vector<vector<string>> up front --
// on multi-million-sample base CSVs (hundreds of columns each) that
// upfront read can need well over 100GB of RSS and gets OOM-killed. Only
// the header/column layout is held for the whole run.
struct DatasetHeader {
    std::vector<std::string> header;
    std::map<std::string, unsigned int> col;
    unsigned int M = 0;
};

static std::vector<std::string> split(const std::string &line, char delim) {
    std::vector<std::string> out;
    std::string cur;
    std::stringstream ss(line);
    while(std::getline(ss, cur, delim))
        out.push_back(cur);
    return out;
}

static std::vector<unsigned int> parseUInts(const std::string &s) {
    std::vector<unsigned int> out;
    for(const std::string &part : split(s, ','))
        if(!part.empty())
            out.push_back((unsigned int)std::stoul(part));
    return out;
}

class StreamingDataset {
public:
    explicit StreamingDataset(const std::string &path) : in_(path.c_str()) {
        if(!in_)
            throw std::runtime_error("cannot open dataset: " + path);
        std::string line;
        if(!std::getline(in_, line))
            throw std::runtime_error("empty dataset: " + path);
        header_.header = split(line, ',');
        for(unsigned int i=0; i<header_.header.size(); i++)
            header_.col[header_.header[i]] = i;
        while(header_.col.count("metric_" + std::to_string(header_.M)))
            header_.M++;
        if(header_.M == 0)
            throw std::runtime_error("dataset has no metric_* columns");
    }

    const DatasetHeader &header() const { return header_; }

    // Reads the next data row into row_; returns false at EOF.
    bool next(std::vector<std::string> &row_out) {
        std::string line;
        if(!std::getline(in_, line))
            return false;
        row_out = split(line, ',');
        return true;
    }

private:
    std::ifstream in_;
    DatasetHeader header_;
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
        else if(a == "--out")
            args.out_csv = need(a);
        else if(a == "--checkpoints")
            args.checkpoints_csv = need(a);
        else if(a == "--max-samples")
            args.max_samples = (unsigned int)std::stoul(need(a));
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if(args.ini_path.empty())
        throw std::runtime_error("--ini is required");
    if(args.dataset_csv.empty())
        throw std::runtime_error("--dataset is required");
    if(args.out_csv.empty())
        throw std::runtime_error("--out is required");
    return args;
}

static double getDouble(const DatasetHeader &h, const std::vector<std::string> &row, const std::string &name) {
    return std::stod(row.at(h.col.at(name)));
}

static int getInt(const DatasetHeader &h, const std::vector<std::string> &row, const std::string &name) {
    return std::stoi(row.at(h.col.at(name)));
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);
        StreamingDataset data(args.dataset_csv);
        const DatasetHeader &h = data.header();
        std::vector<unsigned int> checkpoints = parseUInts(args.checkpoints_csv);

        std::map<std::string, std::map<std::string, std::string> > config = defaultConfig();
        parseConfig(args.ini_path, config);
        config["AdjustPolarDecoder"]["OnlyInit"] = "false";
        config["AdjustPolarDecoder"]["use_AED"] = "true";
        MClassDecoder decoder(config["AdjustPolarDecoder"]);

        const unsigned int n = decoder.codewordLength();
        const unsigned int M = decoder.branchCount();
        if(M != h.M)
            throw std::runtime_error("dataset M does not match decoder branch count");

        std::ofstream out(args.out_csv.c_str());
        if(!out)
            throw std::runtime_error("cannot open output: " + args.out_csv);
        out << std::setprecision(8);
        out << "sample,branch,metric_gap,basin,correct,final_metric";
        for(unsigned int c=0; c<checkpoints.size(); c++) {
            unsigned int step = checkpoints[c];
            out << ",t" << step << "_active"
                << ",t" << step << "_pm_min"
                << ",t" << step << "_pm_gap"
                << ",t" << step << "_pm_mean"
                << ",t" << step << "_pm_max"
                << ",t" << step << "_llr_abs_min"
                << ",t" << step << "_llr_abs_mean"
                << ",t" << step << "_llr_abs_max";
        }
        out << "\n";

        std::vector<std::string> row;
        unsigned int sample_count = 0;
        while((args.max_samples == 0 || sample_count < args.max_samples) && data.next(row)) {
            unsigned int s = sample_count++;
            std::vector<double> y(n, 0.0);
            for(unsigned int i=0; i<n; i++)
                y[i] = getDouble(h, row, "y_" + std::to_string(i));

            std::vector<double> metrics(M, 0.0);
            double best = std::numeric_limits<double>::max();
            for(unsigned int b=0; b<M; b++) {
                metrics[b] = getDouble(h, row, "metric_" + std::to_string(b));
                if(metrics[b] < best)
                    best = metrics[b];
            }

            for(unsigned int b=0; b<M; b++) {
                std::vector<char> candidate;
                double traced_metric = 0.0;
                std::vector<MClassTracePoint> trace;
                bool ok = decoder.decodeBranchTrace(y, b, checkpoints, candidate, traced_metric, trace);
                (void)ok;
                out << s
                    << "," << b
                    << "," << (metrics[b] - best)
                    << "," << ((std::abs(metrics[b] - best) <= 1e-12) ? 1 : 0)
                    << "," << getInt(h, row, "correct_" + std::to_string(b))
                    << "," << metrics[b];
                for(unsigned int c=0; c<checkpoints.size(); c++) {
                    const MClassTracePoint &p = trace.at(c);
                    out << "," << p.active
                        << "," << p.pm_min
                        << "," << p.pm_gap
                        << "," << p.pm_mean
                        << "," << p.pm_max
                        << "," << p.llr_abs_min
                        << "," << p.llr_abs_mean
                        << "," << p.llr_abs_max;
                }
                out << "\n";
            }
        }

        std::cout << "wrote " << args.out_csv << "\n";
        std::cout << "samples=" << sample_count << " M=" << M << " checkpoints=" << checkpoints.size() << "\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
