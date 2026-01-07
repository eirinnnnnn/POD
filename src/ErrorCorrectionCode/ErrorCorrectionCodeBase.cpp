#include "ErrorCorrectionCodeBase.h"
#include "libDebug.h"
#include "libMath.h"
#include "libParser.h"
#include <iostream>
#include <cassert>
#include <typeinfo>

ErrorCorrectionCodeBase::ErrorCorrectionCodeBase(std::map<std::string, std::string> &config){
    // default setting
    ECC_info = "ErrorCorrectionCodeBase";
    encoder_enable = false;
    decoder_enable = false;
    matrix_src = "";
    Gmatrix_path = "";
    Hmatrix_path = "";
    codeword_length = 0;
    message_length = 0;
    dump_Gmatrix_enable = false;
    dump_Hmatrix_enable = false;
    // Grow = 0;

    // deal config
    for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
        if      (pair_idx->first == "matrix_src")    matrix_src = pair_idx->second;
        else if (pair_idx->first == "Gmatrix_path")  Gmatrix_path = pair_idx->second;
        else if (pair_idx->first == "Hmatrix_path")  Hmatrix_path = pair_idx->second;
        else if (pair_idx->first == "dumpGmatrix")   dump_Gmatrix_enable = (pair_idx->second == "True" ? true:false);
        else if (pair_idx->first == "dumpHmatrix")   dump_Hmatrix_enable = (pair_idx->second == "True" ? true:false);

        else WARRING(" key : %s is not find in ErrorCorrectionCodeBase.cpp.",pair_idx->first.c_str());
    }

    // deal Gmatrix and Hmatrix
    // Gmatrix.resize(1);  Gmatrix[0].assign(0,1);
    // Hmatrix.resize(1);  Hmatrix[0].assign(0,1);
    if(matrix_src == "byGmatrix"){
        if(!parseBinaryMatrix(Gmatrix_path,Gmatrix)){
            findNullSpace(Gmatrix, Hmatrix);
        }else{
            WARRING("need to set parity_check.size, message_length, codeword_length");
        }
    }else if(matrix_src == "byHmatrix"){
        if(!parseBinaryMatrix(Hmatrix_path,Hmatrix)){
            std::vector<std::vector<char> > Gmatrix_T, Hmatrix_T;
            transposeMatrix(Hmatrix, Hmatrix_T);
            printMatrix(Hmatrix);
            findNullSpace(Hmatrix_T, Gmatrix_T);
            transposeMatrix(Gmatrix_T, Gmatrix);


        }else{
            WARRING("need to set parity_check.size, message_length, codeword_length");
        }
    }else{
        WARRING("need to set parity_check.size, message_length, codeword_length");
    }
    parity_check.resize(Hmatrix[0].size());
    message_length = Gmatrix.size();
    codeword_length = Gmatrix[0].size();
    // check G * H
    if(DEBUG_MODE){
        std::vector<std::vector<char> > check_matrix;
        matrixMultiplication(Gmatrix,Hmatrix,check_matrix);
        for(unsigned idx=0; idx<check_matrix.size(); idx++)
            for(unsigned jdx=0; jdx<check_matrix[idx].size(); jdx++)
                if(check_matrix[idx][jdx])
                    ERROR("G*H fail");
        INFO("G*H pass");
    }
}

ErrorCorrectionCodeBase::~ErrorCorrectionCodeBase(){}

bool ErrorCorrectionCodeBase::doEncode(std::vector<char> &message, std::vector<char> &codeword){
    matrixMultiplication(message,Gmatrix,codeword);
    return 0;
}

bool ErrorCorrectionCodeBase::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
    assert(received.size() == codeword_length);
    assert(decoded_word.size() == codeword_length);
    for(unsigned int idx=0; idx<codeword_length; idx++)
        decoded_word[idx] = (received[idx] > 0 ? 0:1);
    matrixMultiplication(decoded_word,Hmatrix,parity_check);

    return oneCount(parity_check) == 0;
}

void ErrorCorrectionCodeBase::dumpGmatrix(std::string path){
    char file_name[200];
    sprintf(file_name,"%s/%s_Gmatrix_%d_%d.matrix",path.c_str(),ECC_info.c_str(),message_length,codeword_length);
    FILE *log_data = fopen(file_name,"w");
    printf("%s\n",file_name);
    fprintf(log_data,"%ld %ld\n",Gmatrix.size(),Gmatrix[0].size());
    for(unsigned int idx=0; idx<Gmatrix.size(); idx++){
        for(unsigned int jdx=0; jdx<Gmatrix[0].size(); jdx++)
            fprintf(log_data, "%d",Gmatrix[idx][jdx]);
        fprintf(log_data,"\n");
    }
    fclose(log_data);
}

void ErrorCorrectionCodeBase::dumpHmatrix(std::string path){
    char file_name[200];
    sprintf(file_name,"%s/%s_Hmatrix_%d_%d.matrix",path.c_str(),ECC_info.c_str(),message_length,codeword_length);
    FILE *log_data = fopen(file_name,"w");
    printf("%s\n",file_name);
    fprintf(log_data,"%ld %ld\n",Hmatrix.size(),Hmatrix[0].size());
    for(unsigned int idx=0; idx<Hmatrix.size(); idx++){
        for(unsigned int jdx=0; jdx<Hmatrix[0].size(); jdx++)
            fprintf(log_data, "%d",Hmatrix[idx][jdx]);
        fprintf(log_data,"\n");
    }
    fclose(log_data);
}
