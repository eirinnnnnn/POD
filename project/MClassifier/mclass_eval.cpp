#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Args {
    std::string dataset_csv;
    std::string pred_csv;
    unsigned int topk = 1;
};

struct Dataset {
    std::vector<std::string> header;
    std::map<std::string, unsigned int> col;
    std::vector<std::vector<std::string> > rows;
    unsigned int M = 0;
};

static std::vector<std::string> splitCsvLine(const std::string &line) {
    std::vector<std::string> out;
    std::string cur;
    std::stringstream ss(line);
    while(std::getline(ss, cur, ','))
        out.push_back(cur);
    return out;
}

static Dataset readDataset(const std::string &path) {
    std::ifstream in(path.c_str());
    if(!in)
        throw std::runtime_error("cannot open dataset: " + path);

    Dataset data;
    std::string line;
    if(!std::getline(in, line))
        throw std::runtime_error("empty dataset");
    data.header = splitCsvLine(line);
    for(unsigned int i=0; i<data.header.size(); i++)
        data.col[data.header[i]] = i;

    while(std::getline(in, line)) {
        if(line.empty())
            continue;
        data.rows.push_back(splitCsvLine(line));
    }

    while(data.col.count("metric_" + std::to_string(data.M)))
        data.M++;
    if(data.M == 0)
        throw std::runtime_error("dataset has no metric_* columns");
    return data;
}

static std::vector<std::vector<unsigned int> > readPredictions(const std::string &path) {
    std::ifstream in(path.c_str());
    if(!in)
        throw std::runtime_error("cannot open predictions: " + path);

    std::vector<std::vector<unsigned int> > preds;
    std::string line;
    bool first = true;
    while(std::getline(in, line)) {
        if(line.empty())
            continue;
        std::vector<std::string> fields = splitCsvLine(line);
        if(first && fields.size() && fields[0] == "sample") {
            first = false;
            continue;
        }
        first = false;
        std::vector<unsigned int> row;
        for(unsigned int i=1; i<fields.size(); i++)
            if(!fields[i].empty())
                row.push_back((unsigned int)std::stoul(fields[i]));
        preds.push_back(row);
    }
    return preds;
}

static void usage() {
    std::cout
        << "Usage:\n"
        << "  mclass_eval --dataset dataset.csv --pred predictions.csv [--topk K]\n\n"
        << "Prediction CSV format: sample,idx0,idx1,...\n";
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
        if(a == "--dataset")
            args.dataset_csv = need(a);
        else if(a == "--pred")
            args.pred_csv = need(a);
        else if(a == "--topk")
            args.topk = (unsigned int)std::stoul(need(a));
        else if(a == "-h" || a == "--help") {
            usage();
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + a);
        }
    }
    if(args.dataset_csv.empty())
        throw std::runtime_error("--dataset is required");
    if(args.pred_csv.empty())
        throw std::runtime_error("--pred is required");
    if(args.topk == 0)
        throw std::runtime_error("--topk must be positive");
    return args;
}

static double getDouble(const Dataset &data, const std::vector<std::string> &row, const std::string &name) {
    return std::stod(row.at(data.col.at(name)));
}

static unsigned int getUInt(const Dataset &data, const std::vector<std::string> &row, const std::string &name) {
    return (unsigned int)std::stoul(row.at(data.col.at(name)));
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);
        Dataset data = readDataset(args.dataset_csv);
        std::vector<std::vector<unsigned int> > preds = readPredictions(args.pred_csv);
        if(preds.size() != data.rows.size())
            throw std::runtime_error("prediction row count does not match dataset row count");

        unsigned int teacher_errors = 0;
        unsigned int top1_errors = 0;
        unsigned int topk_errors = 0;
        unsigned int teacher_agree = 0;

        for(unsigned int r=0; r<data.rows.size(); r++) {
            const std::vector<std::string> &row = data.rows[r];
            unsigned int label = getUInt(data, row, "label");
            bool teacher_correct = getUInt(data, row, "teacher_correct") != 0;
            if(!teacher_correct)
                teacher_errors++;

            const std::vector<unsigned int> &p = preds[r];
            if(p.empty())
                throw std::runtime_error("empty prediction row");
            if(p[0] == label)
                teacher_agree++;

            bool top1_correct = false;
            if(p[0] < data.M)
                top1_correct = (getUInt(data, row, "correct_" + std::to_string(p[0])) != 0);
            if(!top1_correct)
                top1_errors++;

            double best_metric = 1e300;
            unsigned int best_idx = p[0];
            unsigned int limit = std::min<unsigned int>(args.topk, p.size());
            for(unsigned int i=0; i<limit; i++) {
                unsigned int idx = p[i];
                if(idx >= data.M)
                    continue;
                double metric = getDouble(data, row, "metric_" + std::to_string(idx));
                if(metric < best_metric) {
                    best_metric = metric;
                    best_idx = idx;
                }
            }
            bool topk_correct = best_idx < data.M &&
                (getUInt(data, row, "correct_" + std::to_string(best_idx)) != 0);
            if(!topk_correct)
                topk_errors++;
        }

        const double N = (double)data.rows.size();
        std::cout << "samples," << data.rows.size() << "\n";
        std::cout << "M," << data.M << "\n";
        std::cout << "teacher_bler," << teacher_errors / N << "\n";
        std::cout << "top1_bler," << top1_errors / N << "\n";
        std::cout << "topk_bler," << topk_errors / N << "\n";
        std::cout << "teacher_agreement," << teacher_agree / N << "\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
