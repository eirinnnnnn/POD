#include "PolarMultiKernal_v2.h"
#include "BCH.h"
#include "libDebug.h"
#include "libParser.h"
#include "libMath.h"
#include <cmath>

PolarMultiKernal_v2::PolarMultiKernal_v2(std::map<std::string, std::string> config) : ErrorCorrectionCodeBase(config){
    // default setting
    ECC_info = "PolarMultiKernal_v2";
    encoder_enable = false;
    decoder_enable = false;
    permutation_src = "";
    dynamic_frozen_process = "frozen";
    // permutation_random_seed = -1;
    kernal_string = "";
    target_raw_BER = 0.1;
    list_size = 0;

    KM = NULL;
    stage = 0;

    receive_order.assign(0,0);
    diverge_flag.assign(0,0);

    // deal config
    for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
        if      (pair_idx->first == "permutation_src")          permutation_src = pair_idx->second;
        else if (pair_idx->first == "dynamic_frozen_process")   dynamic_frozen_process = pair_idx->second;
        else if (pair_idx->first == "target_raw_BER")           target_raw_BER = stod(pair_idx->second);
        else if (pair_idx->first == "list_size")                list_size = stoi(pair_idx->second);
        else if (pair_idx->first == "kernal_string")            kernal_string = pair_idx->second;
        else if (pair_idx->first == "message_length" && message_length == 0)            message_length = stoi(pair_idx->second);
        else {WARRING(" key : %s is not find in PolarMultiKernal_v2.cpp.",pair_idx->first.c_str());}
    }
    // long temp_permutation_random_seed = permutation_random_seed;
    // assert
    if(kernal_string == "")
        ERROR("kernal_string == \"\"");
    KM = new KernalManager(kernal_string);
    stage = KM->get_stage();
    codeword_length = KM->get_codeword_length();
    KM->set_target_BER(target_raw_BER);
    KM->build_Gmatrix(polar_matrix);

    if(DEBUG_MODE){
        printf("--polar_matrix--\n");
        printMatrix(polar_matrix);
        printf("--polar_matrix done-\n");
    }

    if(codeword_length < message_length)
        ERROR("codeword_length(%d) < message_length(%d)",codeword_length,message_length);
    if(list_size < 1)
        ERROR("list_size < 1");

    receive_order.resize(codeword_length);
    Pp.resize(codeword_length);
    PolarGP.resize(codeword_length);
    for(unsigned int row_idx=0; row_idx<codeword_length; row_idx++)
        PolarGP[row_idx].resize(codeword_length);
    diverge_flag.assign(codeword_length,0);
    relation_ship.resize(codeword_length);
    for(unsigned int idx=0; idx<codeword_length; idx++)
        relation_ship[idx].resize(0);

    // get permutation matrix
    permutation_vertor.resize(codeword_length);
    if(permutation_src != ""){
        std::vector<std::vector<char> > permutation_matrix;
        parseBinaryMatrix(permutation_src,permutation_matrix);
        for(unsigned int idx=0; idx<codeword_length; idx++){
            for(unsigned int jdx=0; jdx<codeword_length; jdx++)
                if (permutation_matrix[jdx][idx])
                    permutation_vertor[idx] = jdx;
        }
    }else{
        printf("permutation matrix is idenity\n");
        for(unsigned int idx=0; idx<codeword_length; idx++){
            permutation_vertor[idx] = idx;
        }
    }

    // mem
    list_enable.resize(list_size);
    path_metric.resize(list_size);
    extend_path_metric.resize(list_size*2);
    SCL_mem.resize(list_size);
    for(unsigned int SCL_idx=0; SCL_idx<list_size; SCL_idx++){
        SCL_mem[SCL_idx].resize(stage+1);
        for(unsigned int stage_idx=0; stage_idx<stage+1; stage_idx++)
            SCL_mem[SCL_idx][stage_idx].resize(codeword_length);
    }
    extend_message.resize(codeword_length);
    decoder_enable = true;
    if(DEBUG_MODE){
        INFO("done constructor");
    }
}


PolarMultiKernal_v2::~PolarMultiKernal_v2(){
    delete KM;
}

bool PolarMultiKernal_v2::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
    INFO("doDecode");
    list_enable.assign(list_size,0);
    path_metric.assign(list_size,0);
    // get receive_order, and build Pp
    receive_valid.assign(codeword_length, 1);
    for(unsigned int loop=0; loop<codeword_length; loop++){
        double tmpMin = 99999.9;
        unsigned int Tdx;
        for(unsigned int idx=0; idx<codeword_length; idx++){
            if(receive_valid[idx] && tmpMin > ABS(received[idx])){
                tmpMin = ABS(received[idx]);
                Tdx = idx;
            }
        }
        receive_order[Tdx] = loop;
        receive_valid[Tdx] = 0;
    }
    for(unsigned int idx=0; idx<codeword_length; idx++){
        for(unsigned int jdx=0; jdx<codeword_length; jdx++){
            if (receive_order[jdx] == permutation_vertor[idx]){
                Pp[jdx] = idx;
                break;
            }
        }
    }
    // build relation_ship by Pp
    for(unsigned int idx=0; idx<codeword_length; idx++)
        relation_ship[idx].resize(0);
    for(unsigned int row_idx=0; row_idx<codeword_length; row_idx++)
        for(unsigned int col_idx=0; col_idx<codeword_length; col_idx++)
            PolarGP[row_idx][col_idx] = polar_matrix[row_idx][ Pp[col_idx] ];
    matrixMultiplication(PolarGP, Hmatrix, PolarGPH);
    transposeMatrix(PolarGPH, PolarGPH_T);
    GaussianJordanElimination(PolarGPH_T,1);

    for(unsigned int row_idx=0; row_idx<PolarGPH_T.size(); row_idx++){
        tempV.resize(0);
        for(unsigned int col_idx=0; col_idx<codeword_length; col_idx++){
            if(PolarGPH_T[row_idx][col_idx])
                tempV.push_back(col_idx);
        }

        if(tempV.size() == 0)
            continue;

        relation_ship[tempV[tempV.size()-1]].assign(tempV.begin(), tempV.end());
    }

    // build diverge_flag
    for(unsigned int idx=0; idx<codeword_length; idx++){
        if(relation_ship[idx].size() == 0)   // info
            diverge_flag[idx] = 1;

        if(relation_ship[idx].size()>1 && dynamic_frozen_process == "info")  // dynamic
            diverge_flag[idx] = 1;
    }

    // assign PolarRv, LLR
    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++)
        SCL_mem[0][stage][Pp[codeword_idx]].value = received[codeword_idx]; // need to change to LLR domain
    list_enable[0] = 1;

    for(unsigned int decode_idx=0; decode_idx<codeword_length; decode_idx++){
        INFO("decode_idx %2d",decode_idx);
        // compute node
        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx]){
                KM->do_node_value(decode_idx, 0, SCL_mem[list_idx]);
                INFO("SCL_mem[%2d][0][%2d] = %lf",list_idx,decode_idx,SCL_mem[list_idx][0][decode_idx].value);
            }
        INFO("");

        // diverge or not
        if(diverge_flag[decode_idx])
            infoProcess(decode_idx);
        else
            frozenProcess(decode_idx);

        // update node
        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx])
                KM->update_node_HD(decode_idx, 0, SCL_mem[list_idx]);

        for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
            INFO("path_metric[%2d] = %e", list_idx, path_metric[list_idx]);
        }

        INFO("");
    }

    // check releation ship
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        for(unsigned int codeword_idx=0; codeword_idx<codeword_length && list_enable[list_idx]; codeword_idx++){
            if(diverge_flag[codeword_idx]){
                for(unsigned int relation_idx=0; relation_idx<relation_ship[codeword_idx].size(); relation_idx++)
                    list_enable[list_idx] ^= SCL_mem[list_idx][0][relation_ship[codeword_idx][relation_idx]].HD;
            }
        }
    }

    // pick min path_metric
    double tmp_double = 99999999;
    unsigned int target_list_idx = 0;
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        INFO("path_metric[%2d] = %e", list_idx, path_metric[list_idx]);
        if(list_enable[list_idx] && path_metric[list_idx]<tmp_double){
            target_list_idx = list_idx;
            tmp_double = path_metric[list_idx];
        }
    }
    printf("path_metric,%lf,",path_metric[target_list_idx] );
    // deal output
    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++)
        decoded_word[codeword_idx] = SCL_mem[target_list_idx][stage][Pp[codeword_idx]].HD;

    return true;
}

void PolarMultiKernal_v2::infoProcess(unsigned int decode_idx){
    extend_path_metric.assign(list_size*2,-1);
    unsigned int extend_list_count=0;

    // extend path
    for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
        if(list_enable[list_idx]){
            extend_list_count += 2;
            if(SCL_mem[list_idx][0][decode_idx].value < 0){
                extend_path_metric[list_idx*2  ] = path_metric[list_idx] + ABS(SCL_mem[list_idx][0][decode_idx].value);
                extend_path_metric[list_idx*2+1] = path_metric[list_idx];
            }else{
                extend_path_metric[list_idx*2  ] = path_metric[list_idx];
                extend_path_metric[list_idx*2+1] = path_metric[list_idx] + ABS(SCL_mem[list_idx][0][decode_idx].value);
            }
        }
    }

    // maintain path_metric.size() = list_size, kill large path_metric
    while(extend_list_count > list_size){
        double tmp_double = 0;
        int target_idx = -1;
        for(unsigned int idx=0; idx<extend_path_metric.size(); idx++){
            if(extend_path_metric[idx] > tmp_double){
                tmp_double = extend_path_metric[idx];
                target_idx = idx;
            }
        }
        extend_path_metric[target_idx] = -1;
        extend_list_count--;
    }

    // disable list_enable
    for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
        if(extend_path_metric[list_idx*2] == -1 && extend_path_metric[list_idx*2+1] == -1)
            list_enable[list_idx] = 0;
    }

    // move page
    for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
        if(extend_path_metric[list_idx*2] != -1 && extend_path_metric[list_idx*2+1] != -1){
            for(unsigned int target_idx=0; target_idx<list_enable.size(); target_idx++)
                if(list_enable[target_idx] == 0){
                    clonePage(list_idx,target_idx);
                    list_enable[target_idx] = 1;
                    extend_path_metric[target_idx*2] = extend_path_metric[list_idx*2];
                    extend_path_metric[list_idx*2] = -1;
                    break;
                }
        }
    }

    // set path_metric
    for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
        if(list_enable[list_idx]){
            if(extend_path_metric[list_idx*2] == -1){
                path_metric[list_idx] = extend_path_metric[list_idx*2+1];
                SCL_mem[list_idx][0][decode_idx].HD = 1;
            }else{
                path_metric[list_idx] = extend_path_metric[list_idx*2];
                SCL_mem[list_idx][0][decode_idx].HD = 0;
            }
        }
    }

}

void PolarMultiKernal_v2::frozenProcess(unsigned int decode_idx){
    for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
        if(list_enable[list_idx]){
            SCL_mem[list_idx][0][decode_idx].HD = 0;
            for(unsigned int relation_idx=0; relation_idx<relation_ship[decode_idx].size(); relation_idx++)
                if(decode_idx != relation_ship[decode_idx][relation_idx])
                    SCL_mem[list_idx][0][decode_idx].HD ^= SCL_mem[list_idx][0][relation_ship[decode_idx][relation_idx]].HD;
            if(SIGN(SCL_mem[list_idx][0][decode_idx].value) * (0.5- SCL_mem[list_idx][0][decode_idx].HD) < 0)
                path_metric[list_idx] += ABS(SCL_mem[list_idx][0][decode_idx].value);
        }
    }
}

void PolarMultiKernal_v2::clonePage(unsigned int src, unsigned int dst){
    INFO("clonePage %2d to %2d",src,dst);
    for(unsigned int level_idx=0; level_idx<SCL_mem[src].size(); level_idx++){
        for(unsigned int idx=0; idx<SCL_mem[src][level_idx].size(); idx++){
            SCL_mem[dst][level_idx][idx].value = SCL_mem[src][level_idx][idx].value;
            SCL_mem[dst][level_idx][idx].HD = SCL_mem[src][level_idx][idx].HD;
        }
    }
}
