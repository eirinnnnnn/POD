// Monte-Carlo driver for the exhaustive MLD decoder.
// Structure mirrors project/POD/main.cpp; only the decoder class and the
// config section differ (MLD instead of AdjustPolarDecoder).
#include <iostream>
#include <cstdio>
#include <cstring>
#include "libDebug.h"
#include "libParser.h"
#include "libMath.h"
#include "channelBase.h"
#include "AWGN.h"
#include "ErrorCorrectionCodeBase.h"
#include "MLD.h"

int main(int argc, char const *argv[]){
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

    config["MLD"]["matrix_src"] = "byGmatrix";
    config["MLD"]["Hmatrix_path"] = "";
    config["MLD"]["Gmatrix_path"] = "";

    std::string config_path = "";
    std::string task_folder = "temp_task";
    for(int argc_idx=1; argc_idx<argc; argc_idx++){
        if(!strcmp(argv[argc_idx], "-ini") && argc_idx+1<argc){
            argc_idx++;
            config_path = argv[argc_idx];
            task_folder = config_path.substr(config_path.rfind('/')+1,
                            config_path.rfind('.')-1-config_path.rfind('/'));
            parseConfig(config_path,config);
            continue;
        }
    }
    if(config_path == ""){ printf("usage: -ini <file.ini>\n"); exit(0); }

    int iter_max = std::stoi(config["Monte_Carlo"]["iter_max"]);
    int iter_min = std::stoi(config["Monte_Carlo"]["iter_min"]);
    int error_max = std::stoi(config["Monte_Carlo"]["error_max"]);
    int error_min = std::stoi(config["Monte_Carlo"]["error_min"]);
    int monitor_slot_size = std::stoi(config["Monte_Carlo"]["monitor_slot_size"]);
    long ini_msg_seed = std::stol(config["Monte_Carlo"]["seed_string"]);
    int block_error_count = 0, bit_error_count = 0, current_bit_error = 0;
    int iter_count = 0, noise_bit_count = 0;

    MLD *ECC = new MLD(config["MLD"]);
    channelBase *channel_AWGN = new AWGN(config["AWGN"]);
    channel_AWGN->setCodeRate(ECC->getCodeRate());

    std::vector<char> message(ECC->getMessageLength());
    std::vector<char> codeword(ECC->getCodewordLength());
    std::vector<double> received(ECC->getCodewordLength());
    std::vector<char> decoded_word(ECC->getCodewordLength());
    std::string decode_log;

    system(("mkdir -p "+task_folder).c_str());
    writeBackConfig("./"+task_folder+"/"+task_folder+"_wb.ini",config);

    bool do_monte_carlo = true;
    std::string ini_seed_string = channel_AWGN->getSeedString();
    long msg_seed = ini_msg_seed;
    while(do_monte_carlo){
        iter_count++;
        for(unsigned int i=0;i<ECC->getMessageLength();i++)
            message[i] = (ran0(&msg_seed) > 0.5 ? 1:0);
        ECC->doEncode(message,codeword);
        channel_AWGN->addNoise(codeword,received);
        for(unsigned int idx=0; idx<ECC->getCodewordLength(); idx++){
            if(codeword[idx] && received[idx]>0) noise_bit_count++;
            if(!codeword[idx] && received[idx]<0) noise_bit_count++;
        }
        ECC->doDecode(received,decoded_word,decode_log);
        if(codeword != decoded_word) block_error_count += 1;
        current_bit_error = 0;
        for(unsigned int idx=0; idx<ECC->getCodewordLength(); idx++)
            current_bit_error += codeword[idx]^decoded_word[idx];
        bit_error_count += current_bit_error;

        if(iter_count % monitor_slot_size == 0){
            printf("%s, BLER = %8d/%8d = %8.6f\r",
                channel_AWGN->getChannelInfo().c_str(),
                block_error_count,iter_count,(double)block_error_count/iter_count);
            fflush(stdout);
        }
        if((iter_count >= iter_max && block_error_count >= error_min)
        || (iter_count >= iter_min && block_error_count >= error_max)){
            FILE *log_data = fopen(("./"+task_folder+"/"+"log.txt").c_str(),"a+t");
            printf("%s, BLER = %8d/%8d = %8.6f\n",channel_AWGN->getChannelInfo().c_str(),
                   block_error_count,iter_count,(double)block_error_count/iter_count);
            fflush(stdout);
            fprintf(log_data,"%s, BLER = %8d/%8d = %8.6f\n",
                   channel_AWGN->getChannelInfo().c_str(),
                   block_error_count,iter_count,(double)block_error_count/iter_count);
            fclose(log_data);
            do_monte_carlo = (channel_AWGN->nextChannel() ? true:false);
            channel_AWGN->setSeedString(ini_seed_string);
            msg_seed = ini_msg_seed;
            iter_count = 0; block_error_count = 0; noise_bit_count = 0;
        }
    }
    return 0;
}
