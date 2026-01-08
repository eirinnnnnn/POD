#include <iostream>
#include <cstdio>
#include <cstring>
#include "libDebug.h"
#include "libParser.h"
#include "libMath.h"

#include "channelBase.h"
#include "AWGN.h"

#include "ErrorCorrectionCodeBase.h"
#include "PolarMultiKernal.h"

int main(int argc, char const *argv[]){
    // config default setting
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
    config["PolarMultiKernal"]["Gmatrix_path"] = "../../../data/extend_Golay_12_24_g.matrix";
    config["PolarMultiKernal"]["permutation_src"] = "random";
    config["PolarMultiKernal"]["permutation_random_seed"] = "-1";
    config["PolarMultiKernal"]["dynamic_frozen_process"] = "frozen";
    config["PolarMultiKernal"]["target_raw_BER"] = "0.01";
    config["PolarMultiKernal"]["list_size"] = "4";
    config["PolarMultiKernal"]["kernal_string"] = "23,23,23,753";
    config["PolarMultiKernal"]["bha_value_setting"] = "";
    config["PolarMultiKernal"]["message_length"] = "12";

    config["PolarMultiKernal"]["use_AED"] = "false";
    config["PolarMultiKernal"]["automorphism_src"] = "";
    config["PolarMultiKernal"]["aed_L"] = "0";

    // parse config
    std::string config_path = "";
    std::string task_folder = "temp_task";

auto require_next = [&](int &i){
    if(i + 1 >= argc){
        printf("Missing value for argument: %s\n", argv[i]);
        exit(1);
    }
    i++;
    return std::string(argv[i]);
};

for(int i=1; i<argc; i++){
    std::string arg = argv[i];

    if(arg == "-ini"){
        std::string v = require_next(i);
        config_path = v;

        // default task folder name from ini filename
        {
            size_t slash = config_path.rfind('/');
            size_t dot   = config_path.rfind('.');
            if(dot == std::string::npos) dot = config_path.size();
            if(slash == std::string::npos) slash = (size_t)-1;
            task_folder = config_path.substr(slash+1, dot-(slash+1));
        }

        parseConfig(config_path, config);
        continue;
    }

    if(arg == "-task"){
        task_folder = require_next(i);
        continue;
    }

    if(arg == "-config_example"){
        printf("Generate config_example.ini.\n");
        writeBackConfig("config_example.ini", config);
        exit(0);
    }

    // ===== PolarMultiKernal overrides =====
    if(arg == "-useAED" || arg == "-use_AED"){
        config["PolarMultiKernal"]["useAED"] = require_next(i); // accept 0/1 or true/false
        continue;
    }
    if(arg == "-aed_L"){
        config["PolarMultiKernal"]["aed_L"] = require_next(i);
        continue;
    }
    if(arg == "-automorphism_src"){
        config["PolarMultiKernal"]["automorphism_src"] = require_next(i);
        continue;
    }
    if(arg == "-list_size"){
        config["PolarMultiKernal"]["list_size"] = require_next(i);
        continue;
    }
    if(arg == "-permutation_src"){
        config["PolarMultiKernal"]["permutation_src"] = require_next(i);
        continue;
    }
    if(arg == "-permutation_seed"){
        config["PolarMultiKernal"]["permutation_random_seed"] = require_next(i);
        continue;
    }
    if(arg == "-target_raw_BER"){
        config["PolarMultiKernal"]["target_raw_BER"] = require_next(i);
        continue;
    }
    if(arg == "-kernal_string"){
        config["PolarMultiKernal"]["kernal_string"] = require_next(i);
        continue;
    }

    // ===== AWGN overrides =====
    if(arg == "-snr_start"){
        config["AWGN"]["start"] = require_next(i);
        config["AWGN"]["step_type"] = "SNR";
        continue;
    }
    if(arg == "-snr_end"){
        config["AWGN"]["end"] = require_next(i);
        config["AWGN"]["step_type"] = "SNR";
        continue;
    }
    if(arg == "-snr_step"){
        config["AWGN"]["step"] = require_next(i);
        config["AWGN"]["step_type"] = "SNR";
        continue;
    }
    if(arg == "-awgn_seed"){
        config["AWGN"]["seed_string"] = require_next(i);
        continue;
    }

    // ===== Monte Carlo overrides =====
    if(arg == "-iter_max"){
        config["Monte_Carlo"]["iter_max"] = require_next(i);
        continue;
    }
    if(arg == "-iter_min"){
        config["Monte_Carlo"]["iter_min"] = require_next(i);
        continue;
    }
    if(arg == "-error_max"){
        config["Monte_Carlo"]["error_max"] = require_next(i);
        continue;
    }
    if(arg == "-error_min"){
        config["Monte_Carlo"]["error_min"] = require_next(i);
        continue;
    }
    if(arg == "-monitor"){
        config["Monte_Carlo"]["monitor_slot_size"] = require_next(i);
        continue;
    }
    if(arg == "-msg_seed"){
        config["Monte_Carlo"]["seed_string"] = require_next(i);
        continue;
    }

    // Unknown arg
    printf("Unknown argument: %s\n", argv[i]);
    printf("Hint: use -config_example to generate a template ini.\n");
    exit(1);
}

if(config_path == ""){
    printf("Try to add '-config_example' to get config_example.ini.\n");
    exit(0);
}

    // declare
    int iter_max = std::stoi(config["Monte_Carlo"]["iter_max"]);
    int iter_min = std::stoi(config["Monte_Carlo"]["iter_min"]);
    int error_max = std::stoi(config["Monte_Carlo"]["error_max"]);
    int error_min = std::stoi(config["Monte_Carlo"]["error_min"]);
    int monitor_slot_size = std::stoi(config["Monte_Carlo"]["monitor_slot_size"]);
    long ini_msg_seed = std::stol(config["Monte_Carlo"]["seed_string"]);
    int block_error_count = 0;
    int bit_error_count = 0;
    int current_bit_error = 0;
    int iter_count = 0;
    int noise_bit_count = 0;

    PolarMultiKernal *ECC_PMK = new PolarMultiKernal(config["PolarMultiKernal"]);
    channelBase *channel_AWGN = new AWGN(config["AWGN"]);
    channel_AWGN->setCodeRate(ECC_PMK->getCodeRate());
    // channel_AWGN->setCodeRate(0.5);
    // channel_AWGN->setCodeRate(24.0/47.0);

    std::vector<char> message(ECC_PMK->getMessageLength());
    std::vector<char> codeword(ECC_PMK->getCodewordLength());
    std::vector<double> received(ECC_PMK->getCodewordLength());
    std::vector<char> decoded_word(ECC_PMK->getCodewordLength());
    std::string decode_log;

    // config write back
    system(("mkdir "+task_folder).c_str());
    writeBackConfig("./"+task_folder+"/"+task_folder+"_wb.ini",config);

    // Monte Carlo
    // initialize
    iter_count = 0;
    block_error_count = 0;
    bit_error_count = 0;
    noise_bit_count = 0;
    bool do_monte_carlo = true;
    std::string ini_seed_string = channel_AWGN->getSeedString();
    long msg_seed = ini_msg_seed;
    while(do_monte_carlo){
        iter_count++;
        // random message
        for(unsigned int i_idx=0;i_idx<ECC_PMK->getMessageLength();i_idx++)
            message[i_idx] = (ran0(&msg_seed) > 0.5 ? 1:0);
        // printf("message   ");
        // for(unsigned int i_idx=0; i_idx<message.size(); i_idx++)
        //  printf("%d",message[i_idx]);
        // printf("\n");

        // doEncode
        ECC_PMK->doEncode(message,codeword);
        // printf("codeword  ");
        // for(unsigned int i_idx=0; i_idx<codeword.size(); i_idx++)
        //  printf("%d",codeword[i_idx]);
        // printf("\n");

        // channel
        channel_AWGN->addNoise(codeword,received);
        // received[ECC_PMK->getCodewordLength()-1] = 0.0;
        // for(unsigned int idx=0; idx<ECC_PMK->getCodewordLength(); idx++)
            // printf("(%3d, %d, %6e)\n",idx, codeword[idx], received[idx] );
        for(unsigned int idx=0; idx<ECC_PMK->getCodewordLength(); idx++){
            if(codeword[idx] && received[idx]>0)
                noise_bit_count++;
            if(!codeword[idx] && received[idx]<0)
                noise_bit_count++;
        }

        // printf("received  ");
        // for(unsigned int i_idx=0; i_idx<received.size(); i_idx++)
        //  printf("%d",(received[i_idx] > 0 ?0:1));
        // printf("\n");

        // doDecode
        ECC_PMK->doDecode(received,decoded_word,decode_log);
        // printf("decoded   ");
        // for(unsigned int i_idx=0; i_idx<decoded_word.size(); i_idx++)
        //  printf("%d",decoded_word[i_idx]);
        // printf("\n");

        // deal status
        // printf("ans,%d\n",  codeword != decoded_word ? 1:0);
        if(codeword != decoded_word){
            block_error_count += 1;
            // printf("error in iter_count %d \n",iter_count);
        }

        current_bit_error = 0;
        for(unsigned int idx=0; idx<ECC_PMK->getCodewordLength(); idx++)
            current_bit_error += codeword[idx]^decoded_word[idx];
        bit_error_count += current_bit_error;

        // deal monitor
        if(iter_count % monitor_slot_size == 0){
            printf("%s, BLER = %8d/%8d = %8.6f, BER = %8.6f, RBER = %8.6f\r",
                channel_AWGN->getChannelInfo().c_str(),
                block_error_count,iter_count,
               (double)block_error_count/iter_count,
               (double)bit_error_count/iter_count/ECC_PMK->getCodewordLength(),
               (double)noise_bit_count/iter_count/ECC_PMK->getCodewordLength()
               );
            fflush(stdout);
        }
        // exit(0);

        // deal next loop
        if((iter_count >= iter_max && block_error_count >= error_min) 
        || (iter_count >= iter_min && block_error_count >= error_max)){
            FILE *log_data = fopen(("./"+task_folder+"/"+"log.txt").c_str(),"a+t");
            printf("%s, BLER = %8d/%8d = %8.6f\n",channel_AWGN->getChannelInfo().c_str(),
                                                  block_error_count,iter_count,
                                                  (double)block_error_count/iter_count);
            fflush(stdout);
            fprintf(log_data,"%s, BLER = %8d/%8d = %8.6f\n",channel_AWGN->getChannelInfo().c_str(),
                                                            block_error_count,iter_count,
                                                            (double)block_error_count/iter_count);
            fclose(log_data);
            do_monte_carlo = (channel_AWGN->nextChannel() ? true:false);
            channel_AWGN->setSeedString(ini_seed_string);
            msg_seed = ini_msg_seed;
            iter_count = 0;
            block_error_count = 0;
            noise_bit_count = 0;
        }
    }
    return 0;
}
