#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "libParser.h"
#include "mclass_decoder_wrapper.h"

struct Args {
    std::string ini_path;
    std::string out_csv;
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
        else if(a == "--out")
            args.out_csv = need(a);
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if(args.ini_path.empty())
        throw std::runtime_error("--ini is required");
    if(args.out_csv.empty())
        throw std::runtime_error("--out is required");
    return args;
}

int main(int argc, char **argv) {
    try {
        Args args = parseArgs(argc, argv);
        std::map<std::string, std::map<std::string, std::string> > config = defaultConfig();
        parseConfig(args.ini_path, config);
        config["AdjustPolarDecoder"]["OnlyInit"] = "false";
        config["AdjustPolarDecoder"]["use_AED"] = "true";

        MClassDecoder decoder(config["AdjustPolarDecoder"]);
        std::vector<unsigned int> sizes = decoder.relationSizes();

        std::ofstream out(args.out_csv.c_str());
        if(!out)
            throw std::runtime_error("cannot open output: " + args.out_csv);
        out << "pos,relation_size,bit_type\n";
        for(unsigned int i=0; i<sizes.size(); i++) {
            unsigned int bit_type = sizes[i] > 1 ? 2 : sizes[i];
            out << i << "," << sizes[i] << "," << bit_type << "\n";
        }
        std::cout << "wrote " << args.out_csv << "\n";
        std::cout << "n=" << sizes.size() << "\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "[fatal] " << e.what() << "\n";
        return 1;
    }
}
