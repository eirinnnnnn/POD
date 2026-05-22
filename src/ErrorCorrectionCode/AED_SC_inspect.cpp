// AED_SC_inspect.cpp
// Diagnostic tool for checking whether AED branches genuinely diversify SC/SCL paths.
//
// Purpose:
//   Distinguish
//       no BLER gain
//   from stronger possibilities:
//       same branch-best word,
//       same global AED list union,
//       same valid candidate set,
//       actual received_order collapse.
//
// Usage:
//   ./AED_SC_inspect --config config.ini --samples 10000 --llr-mu 2.0 --llr-sigma 1.0 --seed 1 --out inspect.csv
//
// Notes:
//   - This is not a BLER simulator.
//   - It generates synthetic all-zero-codeword LLRs:
//         L_i ~ Normal(llr_mu, llr_sigma)
//     so positive LLR favors bit 0.
//   - It runs every AED branch on the same received vector.
//   - It records both local branch list size and global union list size across all branches.
//
// Important invariant to check in the output:
//   sample_distinct_all_list_words >= sample_distinct_branch_best_words

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "libMath.h"
#include "libDebug.h"
#include "AED_attemp_20260102.h"

struct BranchCandidate {
    unsigned int branch_id = 0;
    unsigned int list_id = 0;
    bool enabled = false;
    bool valid = false;
    double metric = std::numeric_limits<double>::infinity();
    std::vector<char> word;
};

struct BranchRow {
    unsigned int sample = 0;
    unsigned int branch = 0;

    bool best_enabled = false;
    bool best_valid = false;
    double best_metric = std::numeric_limits<double>::infinity();
    uint64_t best_hash = 0;
    std::string best_word;

    unsigned int local_enabled_count = 0;
    unsigned int local_valid_count = 0;
    unsigned int local_distinct_list_words = 0;
    unsigned int local_distinct_valid_list_words = 0;
};

static std::string trim(const std::string &s){
    size_t a = 0;
    while(a < s.size() && std::isspace((unsigned char)s[a])) a++;

    size_t b = s.size();
    while(b > a && std::isspace((unsigned char)s[b - 1])) b--;

    return s.substr(a, b - a);
}

static std::map<std::string, std::string> read_ini_simple(const std::string &path){
    std::ifstream fin(path.c_str());
    if(!fin)
        throw std::runtime_error("cannot open config: " + path);

    std::map<std::string, std::string> cfg;
    std::string line;
    unsigned int line_no = 0;

    while(std::getline(fin, line)){
        line_no++;

        // Simple comment stripping.
        // Avoid using # or ; inside values.
        size_t sharp = line.find('#');
        if(sharp != std::string::npos)
            line = line.substr(0, sharp);

        size_t semi = line.find(';');
        if(semi != std::string::npos)
            line = line.substr(0, semi);

        line = trim(line);
        if(line.empty())
            continue;

        // Ignore section headers.
        if(line.front() == '[' && line.back() == ']')
            continue;

        size_t eq = line.find('=');
        if(eq == std::string::npos){
            std::cerr << "[warn] ignore config line " << line_no << ": " << line << "\n";
            continue;
        }

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));
        cfg[key] = val;
    }

    return cfg;
}

static std::string bits_to_string(const std::vector<char> &v){
    std::string s;
    s.reserve(v.size());
    for(char b : v)
        s.push_back(b ? '1' : '0');
    return s;
}

static std::string uint_vec_to_string(const std::vector<unsigned int> &v){
    std::ostringstream oss;
    oss << "[";
    for(size_t i = 0; i < v.size(); i++){
        if(i) oss << ",";
        oss << v[i];
    }
    oss << "]";
    return oss.str();
}

static std::string relation_to_string(const std::vector<unsigned int> &v){
    if(v.empty())
        return "info";

    std::ostringstream oss;
    oss << "{";
    for(size_t i = 0; i < v.size(); i++){
        if(i) oss << ",";
        oss << v[i];
    }
    oss << "}";
    return oss.str();
}

static uint64_t fnv1a_hash_uvec(const std::vector<unsigned int> &v){
    uint64_t h = 1469598103934665603ULL;
    for(unsigned int x : v){
        h ^= (uint64_t)x + 0x9e3779b97f4a7c15ULL;
        h *= 1099511628211ULL;
    }
    return h;
}

static uint64_t fnv1a_hash_word(const std::vector<char> &v){
    uint64_t h = 1469598103934665603ULL;
    for(char x : v){
        h ^= (uint64_t)(unsigned char)x;
        h *= 1099511628211ULL;
    }
    return h;
}

class AEDInspector : public AdjustPolarDecoder {
public:
    explicit AEDInspector(std::map<std::string, std::string> config)
        : AdjustPolarDecoder(config) {}

    unsigned int n() const {
        return codeword_length;
    }

    unsigned int L() const {
        return list_size;
    }

    unsigned int effective_branch_count() const {
        if(received_order_set.empty())
            return 1;
        return (unsigned int)received_order_set.size();
    }

    const std::vector<unsigned int>& effective_branch_order(unsigned int i) const {
        if(received_order_set.empty()){
            if(i != 0)
                throw std::runtime_error("branch index out of range");
            return received_order;
        }
        return received_order_set.at(i);
    }

    void dump_structure(std::ostream &os){
        os << "==== decoder structure ====\n";
        os << "n=" << codeword_length
           << " k=" << message_length
           << " list_size=" << list_size
           << " branches=" << effective_branch_count() << "\n";

        os << "\n-- received_order base --\n";
        os << uint_vec_to_string(received_order) << "\n";

        os << "\n-- branch received_order hashes --\n";
        std::set<uint64_t> order_hashes;

        for(unsigned int i = 0; i < effective_branch_count(); i++){
            const std::vector<unsigned int> &ord = effective_branch_order(i);
            uint64_t h = fnv1a_hash_uvec(ord);
            order_hashes.insert(h);

            os << "branch " << i
               << " hash=" << h
               << " order=" << uint_vec_to_string(ord) << "\n";
        }

        os << "distinct_received_orders=" << order_hashes.size() << "\n";

        os << "\n-- relation_ship / causal frozen relations --\n";
        for(unsigned int i = 0; i < codeword_length; i++){
            os << std::setw(3) << i
               << " : " << relation_to_string(relation_ship[i])
               << " diverge=" << (int)diverge_flag[i] << "\n";
        }

        os << "==== end structure ====\n";
    }

    bool decode_order_collect_all(const std::vector<double> &received,
                                  const std::vector<unsigned int> &order,
                                  unsigned int branch_id,
                                  std::vector<BranchCandidate> &out){
        out.clear();

        if(order.size() != codeword_length){
            std::cerr << "[warn] order size mismatch\n";
            return false;
        }

        // Reset working memory for all list entries.
        // This mirrors decode_with_order_scl().
        for(unsigned int list_idx = 0; list_idx < list_size; list_idx++){
            for(unsigned int level_idx = 0; level_idx < SCL_mem[list_idx].size(); level_idx++){
                for(unsigned int idx = 0; idx < SCL_mem[list_idx][level_idx].size(); idx++){
                    SCL_mem[list_idx][level_idx][idx].HD = 0;
                    SCL_mem[list_idx][level_idx][idx].value = 0.0;
                }
            }
        }

        std::fill(list_enable.begin(), list_enable.end(), 0);
        std::fill(path_metric.begin(), path_metric.end(), 0.0);
        list_enable[0] = 1;

        // Apply received LLRs through the branch order.
        for(unsigned int codeword_idx = 0; codeword_idx < received.size(); codeword_idx++){
            unsigned int target = order[codeword_idx];

            if(target >= codeword_length){
                std::cerr << "[warn] target index out of range: " << target << "\n";
                return false;
            }

            if(bha_value_setting[codeword_idx] == '?'){
                SCL_mem[0][stage][target].value = received[codeword_idx];
            }else if(bha_value_setting[codeword_idx] == '1'){
                SCL_mem[0][stage][target].value = 99999999999999999.9999;
            }else if(bha_value_setting[codeword_idx] == '0'){
                SCL_mem[0][stage][target].value = 0.0;
            }else{
                throw std::runtime_error("bad bha_value_setting");
            }
        }

        // SCL recursion, following the existing decoder.
        for(unsigned int decode_idx = 0; decode_idx < codeword_length; decode_idx++){
            for(unsigned int list_idx = 0; list_idx < list_size; list_idx++){
                if(list_enable[list_idx])
                    do_node_value(decode_idx, 0, SCL_mem[list_idx]);
            }

            if(diverge_flag[decode_idx])
                infoProcess(decode_idx);
            else
                frozenProcess(decode_idx);

            for(unsigned int list_idx = 0; list_idx < list_size; list_idx++){
                if(list_enable[list_idx])
                    update_node_HD(decode_idx, 0, SCL_mem[list_idx]);
            }
        }

        // Same relation check as decode_with_order_scl().
        for(unsigned int list_idx = 0; list_idx < list_size; list_idx++){
            for(unsigned int codeword_idx = 0;
                codeword_idx < codeword_length && list_enable[list_idx];
                codeword_idx++){

                if(diverge_flag[codeword_idx]){
                    for(unsigned int relation_idx = 0;
                        relation_idx < relation_ship[codeword_idx].size();
                        relation_idx++){

                        list_enable[list_idx] ^= SCL_mem[list_idx][0][relation_ship[codeword_idx][relation_idx]].HD;
                    }
                }
            }
        }

        // Collect every enabled SCL path from this branch.
        for(unsigned int list_idx = 0; list_idx < list_size; list_idx++){
            BranchCandidate cand;
            cand.branch_id = branch_id;
            cand.list_id = list_idx;
            cand.enabled = (list_enable[list_idx] != 0);
            cand.metric = path_metric[list_idx];
            cand.word.assign(codeword_length, 0);

            if(cand.enabled){
                for(unsigned int codeword_idx = 0; codeword_idx < received.size(); codeword_idx++){
                    cand.word[codeword_idx] =
                        SCL_mem[list_idx][stage][order[codeword_idx]].HD;
                }

                std::vector<char> syndrome;
                matrixMultiplication(cand.word, Hmatrix, syndrome);
                cand.valid = (oneCount(syndrome) == 0);
            }

            out.push_back(cand);
        }

        return true;
    }
};

struct Args {
    std::string config;
    std::string out_csv;
    unsigned int samples = 1000;
    unsigned int seed = 1;
    double llr_mu = 2.0;
    double llr_sigma = 1.0;
    unsigned int max_print_samples = 5;
};

static Args parse_args(int argc, char **argv){
    Args args;

    for(int i = 1; i < argc; i++){
        std::string a = argv[i];

        auto need = [&](const std::string &name)->std::string{
            if(i + 1 >= argc)
                throw std::runtime_error("missing value after " + name);
            return argv[++i];
        };

        if(a == "--config"){
            args.config = need(a);
        }else if(a == "--out"){
            args.out_csv = need(a);
        }else if(a == "--samples"){
            args.samples = (unsigned int)std::stoul(need(a));
        }else if(a == "--seed"){
            args.seed = (unsigned int)std::stoul(need(a));
        }else if(a == "--llr-mu"){
            args.llr_mu = std::stod(need(a));
        }else if(a == "--llr-sigma"){
            args.llr_sigma = std::stod(need(a));
        }else if(a == "--max-print-samples"){
            args.max_print_samples = (unsigned int)std::stoul(need(a));
        }else if(a == "--help" || a == "-h"){
            std::cout
                << "Usage:\n"
                << "  ./AED_SC_inspect --config config.ini "
                << "[--samples 1000] [--llr-mu 2.0] [--llr-sigma 1.0] "
                << "[--seed 1] [--out inspect.csv]\n";
            std::exit(0);
        }else{
            throw std::runtime_error("unknown arg: " + a);
        }
    }

    if(args.config.empty())
        throw std::runtime_error("--config is required");

    return args;
}

int main(int argc, char **argv){
    try{
        Args args = parse_args(argc, argv);
        std::map<std::string, std::string> cfg = read_ini_simple(args.config);

        // Force AED on if automorphism_src is supplied but use_AED is missing.
        // Comment this out if you want strict config behavior.
        if(cfg.count("automorphism_src") && !cfg.count("use_AED"))
            cfg["use_AED"] = "true";

        AEDInspector dec(cfg);
        dec.dump_structure(std::cout);

        std::ofstream fout;
        if(!args.out_csv.empty()){
            fout.open(args.out_csv.c_str());
            if(!fout)
                throw std::runtime_error("cannot open output csv: " + args.out_csv);

            fout
                << "sample,branch,"
                << "best_enabled,best_valid,best_metric,best_hash,best_word,"
                << "local_enabled_count,local_valid_count,"
                << "local_distinct_list_words,local_distinct_valid_list_words,"
                << "sample_enabled_count,sample_valid_count,"
                << "sample_distinct_branch_best_words,"
                << "sample_distinct_all_list_words,"
                << "sample_distinct_valid_list_words,"
                << "sample_branch_best_all_same,"
                << "sample_best_branch_differs_from_branch0,"
                << "sample_list_union_size_leq_L\n";
        }

        std::mt19937_64 rng(args.seed);
        std::normal_distribution<double> gauss(args.llr_mu, args.llr_sigma);

        uint64_t total_enabled_lists = 0;
        uint64_t total_valid_lists = 0;
        uint64_t total_distinct_branch_best = 0;
        uint64_t total_distinct_all_lists = 0;
        uint64_t total_distinct_valid_lists = 0;
        uint64_t samples_best_branch_differs_from_0 = 0;
        uint64_t samples_branch_best_all_same = 0;
        uint64_t samples_list_union_size_leq_L = 0;

        for(unsigned int s = 0; s < args.samples; s++){
            std::vector<double> received(dec.n(), 0.0);
            for(double &x : received)
                x = gauss(rng);

            std::vector<BranchRow> rows;
            rows.reserve(dec.effective_branch_count());

            std::set<std::string> sample_branch_best_words;
            std::set<std::string> sample_all_list_words;
            std::set<std::string> sample_valid_list_words;

            unsigned int sample_enabled_count = 0;
            unsigned int sample_valid_count = 0;

            std::vector<BranchCandidate> branch_best(dec.effective_branch_count());
            std::vector<char> branch_has_best(dec.effective_branch_count(), 0);

            for(unsigned int b = 0; b < dec.effective_branch_count(); b++){
                std::vector<BranchCandidate> list_cands;
                bool ok = dec.decode_order_collect_all(
                    received,
                    dec.effective_branch_order(b),
                    b,
                    list_cands
                );

                BranchRow row;
                row.sample = s;
                row.branch = b;

                if(!ok){
                    rows.push_back(row);
                    continue;
                }

                BranchCandidate best;
                best.branch_id = b;
                best.metric = std::numeric_limits<double>::infinity();
                bool have_best = false;

                std::set<std::string> local_words;
                std::set<std::string> local_valid_words;

                for(const BranchCandidate &c : list_cands){
                    if(!c.enabled)
                        continue;

                    row.local_enabled_count++;
                    sample_enabled_count++;

                    std::string w = bits_to_string(c.word);

                    local_words.insert(w);
                    sample_all_list_words.insert(w);

                    if(c.valid){
                        row.local_valid_count++;
                        sample_valid_count++;

                        local_valid_words.insert(w);
                        sample_valid_list_words.insert(w);
                    }

                    if(!have_best || c.metric < best.metric){
                        best = c;
                        have_best = true;
                    }
                }

                row.local_distinct_list_words = (unsigned int)local_words.size();
                row.local_distinct_valid_list_words = (unsigned int)local_valid_words.size();

                if(have_best){
                    branch_best[b] = best;
                    branch_has_best[b] = 1;

                    row.best_enabled = true;
                    row.best_valid = best.valid;
                    row.best_metric = best.metric;
                    row.best_hash = fnv1a_hash_word(best.word);
                    row.best_word = bits_to_string(best.word);

                    sample_branch_best_words.insert(row.best_word);
                }

                rows.push_back(row);
            }

            bool sample_branch_best_all_same =
                (sample_branch_best_words.size() == 1);

            bool sample_best_branch_differs_from_branch0 = false;
            if(dec.effective_branch_count() > 1 && branch_has_best[0]){
                std::string w0 = bits_to_string(branch_best[0].word);

                for(unsigned int b = 1; b < dec.effective_branch_count(); b++){
                    if(branch_has_best[b] &&
                       bits_to_string(branch_best[b].word) != w0){

                        sample_best_branch_differs_from_branch0 = true;
                        break;
                    }
                }
            }

            bool sample_list_union_size_leq_L =
                (sample_all_list_words.size() <= dec.L());

            total_enabled_lists += sample_enabled_count;
            total_valid_lists += sample_valid_count;
            total_distinct_branch_best += sample_branch_best_words.size();
            total_distinct_all_lists += sample_all_list_words.size();
            total_distinct_valid_lists += sample_valid_list_words.size();

            if(sample_branch_best_all_same)
                samples_branch_best_all_same++;

            if(sample_best_branch_differs_from_branch0)
                samples_best_branch_differs_from_0++;

            if(sample_list_union_size_leq_L)
                samples_list_union_size_leq_L++;

            // Now write CSV rows after sample-level global union sizes are known.
            if(fout){
                for(const BranchRow &row : rows){
                    fout
                        << row.sample << ","
                        << row.branch << ","
                        << (row.best_enabled ? 1 : 0) << ","
                        << (row.best_valid ? 1 : 0) << ",";

                    if(row.best_enabled){
                        fout
                            << std::setprecision(17) << row.best_metric << ","
                            << row.best_hash << ","
                            << row.best_word << ",";
                    }else{
                        fout << "inf,0,,";
                    }

                    fout
                        << row.local_enabled_count << ","
                        << row.local_valid_count << ","
                        << row.local_distinct_list_words << ","
                        << row.local_distinct_valid_list_words << ","
                        << sample_enabled_count << ","
                        << sample_valid_count << ","
                        << sample_branch_best_words.size() << ","
                        << sample_all_list_words.size() << ","
                        << sample_valid_list_words.size() << ","
                        << (sample_branch_best_all_same ? 1 : 0) << ","
                        << (sample_best_branch_differs_from_branch0 ? 1 : 0) << ","
                        << (sample_list_union_size_leq_L ? 1 : 0)
                        << "\n";
                }
            }

            if(s < args.max_print_samples){
                std::cout << "\n[sample " << s << "]\n";
                std::cout
                    << "sample_distinct_branch_best_words="
                    << sample_branch_best_words.size()
                    << " sample_distinct_all_list_words="
                    << sample_all_list_words.size()
                    << " sample_distinct_valid_list_words="
                    << sample_valid_list_words.size()
                    << "\n";

                for(unsigned int b = 0; b < dec.effective_branch_count(); b++){
                    std::cout << "  branch " << b;

                    if(branch_has_best[b]){
                        std::cout
                            << " best_metric=" << branch_best[b].metric
                            << " valid=" << branch_best[b].valid
                            << " word=" << bits_to_string(branch_best[b].word)
                            << "\n";
                    }else{
                        std::cout << " no_enabled_path\n";
                    }
                }
            }
        }

        double S = (double)args.samples;

        std::cout << "\n==== aggregate diversity report ====\n";
        std::cout
            << "samples=" << args.samples
            << " branches=" << dec.effective_branch_count()
            << " list_size=" << dec.L()
            << "\n";

        std::cout
            << "avg_enabled_lists_per_sample="
            << (total_enabled_lists / S)
            << "\n";

        std::cout
            << "avg_valid_lists_per_sample="
            << (total_valid_lists / S)
            << "\n";

        std::cout
            << "avg_distinct_branch_best_words="
            << (total_distinct_branch_best / S)
            << "\n";

        std::cout
            << "avg_distinct_all_list_words="
            << (total_distinct_all_lists / S)
            << "\n";

        std::cout
            << "avg_distinct_valid_list_words="
            << (total_distinct_valid_lists / S)
            << "\n";

        std::cout
            << "Pr(branch_best_all_same)="
            << (samples_branch_best_all_same / S)
            << "\n";

        std::cout
            << "Pr(best_branch_differs_from_branch0)="
            << (samples_best_branch_differs_from_0 / S)
            << "\n";

        std::cout
            << "Pr(list_union_size <= L)="
            << (samples_list_union_size_leq_L / S)
            << "\n";

        std::cout << "==== end report ====\n";

        if(fout)
            std::cout << "[csv saved] " << args.out_csv << "\n";

    }catch(const std::exception &e){
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }

    return 0;
}