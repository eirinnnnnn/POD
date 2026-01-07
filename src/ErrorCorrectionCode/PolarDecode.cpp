#include "PolarDecode.h"
#include "libMath.h"
#include "libDebug.h"
#include "libParser.h"
#include "BCH.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>


PolarDecode::PolarDecode(std::map<std::string, std::string> config) : PolarCode(config){
    // default setting
    ECC_info = "PolarDecode";
    encoder_enable = false;
    decoder_enable = false;

    permutation_src = "";
    dynamic_frozen_process = "frozen";
    permutation_random_seed = -1;
    // deal config
    for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
        if      (pair_idx->first == "permutation_src")          permutation_src = pair_idx->second;
        else if (pair_idx->first == "dynamic_frozen_process")   dynamic_frozen_process = pair_idx->second;
        else if (pair_idx->first == "permutation_random_seed")  permutation_random_seed = stol(pair_idx->second);
        else WARRING(" key : %s is not find in PolarDecode.cpp.",pair_idx->first.c_str());
    }

    permutation_array.resize(codeword_length);
    // internal_received.resize(codeword_length);
    // internal_decoded_word.resize(codeword_length);

    if(matrix_src == "byFile"){
        unsigned int Prow = codeword_length;
        unsigned int Pcol = codeword_length;

        if(permutation_src == "random")
            permutationMatrix(permutation_matrix,codeword_length,&permutation_random_seed);
        else if(permutation_src != ""){
            parseBinaryMatrix(permutation_src,permutation_matrix);
            Prow = permutation_matrix.size();
            Pcol = permutation_matrix[0].size();
        }else{
            permutation_matrix.resize(codeword_length);
            for(unsigned int idx=0; idx<codeword_length; idx++){
                permutation_matrix[idx].assign(codeword_length,0);
                permutation_matrix[idx][idx] = 1;
            }
        }
        if(Pcol != codeword_length || Prow != codeword_length){
            ERROR("Pcol != codeword_length || Prow != codeword_length");    exit(1);
        }

        for(unsigned int idx=0; idx<codeword_length; idx++)
            for(unsigned int jdx=0; jdx<codeword_length; jdx++)
                if(permutation_matrix[idx][jdx])
                    permutation_array[jdx] = idx;
        buildRelationShip();
        decoder_enable = true;
    }else if(matrix_src == "byBCH"){
        std::map<std::string, std::string> BCH_config;
        BCH_config["matrix_src"] = "byAlgorithm";
        BCH_config["codeword_length"] = std::to_string(codeword_length-1);
        BCH_config["message_length"] = std::to_string(message_length);
        BCH ECC_BCH = BCH(BCH_config);
        std::vector<std::vector<char> >BCH_Hmatrix = ECC_BCH.getHmatrix();
        Hmatrix.resize(BCH_Hmatrix.size()+1);
        for(unsigned int idx=0; idx<BCH_Hmatrix.size(); idx++){
            Hmatrix[idx].assign(BCH_Hmatrix[idx].size()+1,0);
            for(unsigned int jdx=0; jdx<BCH_Hmatrix[idx].size(); jdx++)
                Hmatrix[idx][jdx] = BCH_Hmatrix[idx][jdx];
            Hmatrix[idx][BCH_Hmatrix[idx].size()] = 1;
        }
        Hmatrix[BCH_Hmatrix.size()].assign(BCH_Hmatrix[0].size()+1,0);
        Hmatrix[BCH_Hmatrix.size()][BCH_Hmatrix[0].size()] = 1;
        
        std::vector<unsigned int> BCH_DtoA = ECC_BCH.getDtoA();
        for(unsigned int idx=0; idx<codeword_length-1; idx++)
            permutation_array[idx] = codeword_length -1 - BCH_DtoA[idx];
        permutation_array[codeword_length-1] = codeword_length-1;
        permutationMatrix(permutation_matrix,permutation_array);
        
        // printMatrix(Hmatrix);
        // printMatrix(permutation_matrix);

        buildRelationShip();
        decoder_enable = true;
    }


    // printf("message_length %d\n", message_length);
    // printf("codeword_length %d\n", codeword_length);
    // printf("decoder_enable %d\n", decoder_enable);
    // printf("encoder_enable %d\n", encoder_enable);
    // printf("done\n");
    // exit(1);
}

PolarDecode::~PolarDecode(){}

bool PolarDecode::buildRelationShip(){
    std::vector<std::vector<char> > temp_matrixA,temp_matrixB;
    std::vector<std::vector<char> > valid_matrix;
    std::vector<std::vector<char> > Polar_code_matrix;
    Polar_code_matrix.resize(codeword_length);
    for(unsigned int idx=0; idx < codeword_length; idx++){
        Polar_code_matrix[idx].assign(codeword_length,0);
        for(unsigned int sub_idx=0; sub_idx < codeword_length; sub_idx++)
            if((sub_idx | idx) == idx)
                Polar_code_matrix[idx][sub_idx] = 1;
    }   

    matrixMultiplication(permutation_matrix,Hmatrix,temp_matrixA);
    matrixMultiplication(Polar_code_matrix,temp_matrixA,temp_matrixB);
    transposeMatrix(temp_matrixB,valid_matrix);
    GaussianJordanElimination(valid_matrix,1);
    // printMatrix(valid_matrix);
    for(unsigned int idx=0; idx<codeword_length; idx++)
        relation_ship[idx].resize(0);
    for(unsigned int idx=0; idx<valid_matrix.size(); idx++){
        int target_idx = -1;
        for(unsigned int jdx=0; jdx<codeword_length; jdx++){
            if(valid_matrix[idx][codeword_length-1-jdx]){
                if(target_idx == -1)
                    target_idx = codeword_length-1-jdx;
                relation_ship[target_idx].push_back(codeword_length-1-jdx);
            }
        }
    }
    for(unsigned int idx=0; idx<codeword_length; idx++){
        received_order[idx] = permutation_array[idx];
        if(relation_ship[idx].size() == 0)
            diverge_flag[idx] = 1;
        else if(relation_ship[idx].size() > 1 && dynamic_frozen_process == "info")
            diverge_flag[idx] = 1;
        else
            diverge_flag[idx] = 0;
        // printf("%2d %d %d\n",idx,diverge_flag[idx],relation_ship[idx].size() );
    }

    return 0;
}
