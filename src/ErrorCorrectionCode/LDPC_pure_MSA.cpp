#include "LDPC_pure_MSA.h"
#include "libDebug.h"
#include "libMath.h"
#include <cmath>
#include <cassert>
#include <cstdio>

void PureMSA::printf_BN(){
    printf("BN INFO\n");
    for (unsigned int idx=0; idx<BN.size(); idx++){
        printf("%3d ",idx );
        for (unsigned int jdx=0; jdx<BN[idx].link.size(); jdx++)
            printf("(%3d, %+8.3f)", BN[idx].link[jdx], BN[idx].Q[jdx]);
        printf("\n");
    }
}

void PureMSA::printf_CN(){
    printf("CN INFO\n");
    for (unsigned int idx=0; idx<CN.size(); idx++){
        printf("%3d %d ",idx , CN[idx].allSign);
        for (unsigned int jdx=0; jdx<CN[idx].link.size(); jdx++)
            printf("(%3d, %+8.3f)", CN[idx].link[jdx], CN[idx].R[jdx]);
        printf("\n");
    }
}

PureMSA::PureMSA(std::map<std::string, std::string> &config) : ErrorCorrectionCodeBase(config){
    // default setting
    Maxiter = 0;
    // deal config
    for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
        if      (pair_idx->first == "Maxiter")                  Maxiter = stoi(pair_idx->second);
        // else if (pair_idx->first == "target_raw_BER")           target_raw_BER = stod(pair_idx->second);
        else {WARRING(" key : %s is not find in PolarMultiKernal.cpp.",pair_idx->first.c_str());}
    }

    // assert 
    if(matrix_src != "byHmatrix")
        ERROR("LDPC should initail with loading H matrix");

    // set BN, CN
    BN.resize(Hmatrix.size());
    CN.resize(Hmatrix[0].size());
    for(unsigned int idx=0; idx<Hmatrix.size(); idx++)
        for(unsigned int jdx=0; jdx<Hmatrix[idx].size(); jdx++){
            if(Hmatrix[idx][jdx]){
                BN[idx].link.push_back(jdx);
                CN[jdx].link.push_back(idx);
            }
        }
    for(unsigned int idx=0; idx<BN.size(); idx++){
        BN[idx].Q.resize(BN[idx].link.size());
    }    
    for(unsigned int idx=0; idx<CN.size(); idx++){
        CN[idx].R.resize(CN[idx].link.size());
    }

    if(DEBUG_MODE){
        printf("link BN -> CN\n");
        for(unsigned int idx=0; idx<BN.size(); idx++){
            printf("BN %3d -> ", idx);
            for(unsigned int jdx=0; jdx<BN[idx].link.size(); jdx++){
                printf("%3d, ",BN[idx].link[jdx] );
            }
            printf("\n");
        }

        printf("link CN -> BN\n");
        for(unsigned int idx=0; idx<CN.size(); idx++){
            printf("CN %3d -> ", idx);
            for(unsigned int jdx=0; jdx<CN[idx].link.size(); jdx++){
                printf("%3d, ",CN[idx].link[jdx] );
            }
            printf("\n");
        }

    }

    printf("Maxiter %d\n",Maxiter );
}

PureMSA::~PureMSA(){}

bool PureMSA::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
    
    // initial 
    for (unsigned int idx=0; idx<BN.size(); idx++){
        BN[idx].received = received[idx];
        for (unsigned int jdx=0; jdx<BN[idx].Q.size(); jdx++)
            BN[idx].Q[jdx] = BN[idx].received;
    }
    for(unsigned int idx=0; idx<CN.size(); idx++){
        for (unsigned int jdx=0; jdx<CN[idx].R.size(); jdx++) 
            CN[idx].R[jdx] = 0;
    }

    if(DEBUG_MODE){
        printf_BN();
        printf_CN();
    }

    for(unsigned int iter=0; iter<Maxiter; iter++){
        
        CNP_SPA();
        if(DEBUG_MODE){
            printf_CN(); 
        }
               
        BNP();
        if(DEBUG_MODE){
            printf_BN();
        }

        if (ParityCheck(decoded_word)){
            break;
        }
    }
    return 0;
}

void PureMSA::BNP(){
    for (unsigned int idx=0; idx<BN.size(); idx++){
        double allR_L = BN[idx].received;
        for (unsigned int jdx=0; jdx<BN[idx].link.size(); jdx++){
            checkNode *linkedCN = &(CN[BN[idx].link[jdx]]);
            for (unsigned int cidx=0; cidx<linkedCN->link.size(); cidx++)
                if(linkedCN->link[cidx] == idx){
                    BN[idx].Q[jdx] = linkedCN->R[cidx];
                    allR_L += BN[idx].Q[jdx];
                    break;
                }
        }

        for (unsigned int jdx=0; jdx<BN[idx].Q.size(); jdx++){
            BN[idx].Q[jdx] = allR_L - BN[idx].Q[jdx];
        }
    }
}

char PureMSA::CNP_MSA(){
    char ET = 1;
    for (unsigned int idx=0; idx<CN.size(); idx++){
        double min = 999999999;
        double secMin = 999999999;
        unsigned int idx_min=0;
        CN[idx].allSign = 1;
        for (unsigned int jdx=0; jdx<CN[idx].link.size(); jdx++){
            bitNode *linkedBN = &(BN[CN[idx].link[jdx]]);
            for (unsigned int bidx=0; bidx<linkedBN->link.size(); bidx++)
                if(linkedBN->link[bidx] == idx){
                    double Q = linkedBN->Q[bidx];
                    CN[idx].R[jdx] = SIGN(Q);
                    Q = ABS(Q);
                    CN[idx].allSign *= CN[idx].R[jdx];
                    
                    if (Q < min){
                        secMin = min;
                        min = Q;
                        idx_min = jdx;
                    }else if (Q < secMin){
                        secMin = Q;
                    }
                    break;
                }
        }
        ET &= (CN[idx].allSign == 1);
        for (unsigned int jdx=0; jdx<CN[idx].R.size(); jdx++){
            CN[idx].R[jdx] *= CN[idx].allSign * (jdx == idx_min ? secMin:min);
        }
    }

    return ET == 1;
}

char PureMSA::CNP_SPA(){
    char ET = 1;
    for (unsigned int idx=0; idx<CN.size(); idx++){
        double all_SPA_value = 1;
        CN[idx].allSign = 1;
        for (unsigned int jdx=0; jdx<CN[idx].link.size(); jdx++){
            bitNode *linkedBN = &(BN[CN[idx].link[jdx]]);
            for (unsigned int bidx=0; bidx<linkedBN->link.size(); bidx++)
                if(linkedBN->link[bidx] == idx){
                    double Q = linkedBN->Q[bidx];
                    CN[idx].R[jdx] = tanh(Q/2.0);
                    all_SPA_value *= CN[idx].R[jdx];
                    CN[idx].allSign *= SIGN(Q);
                    break;
                }
        }
        ET &= (CN[idx].allSign == 1);
        for (unsigned int jdx=0; jdx<CN[idx].R.size(); jdx++){
            CN[idx].R[jdx] = 2*atanh(all_SPA_value / CN[idx].R[jdx]);
        }
    }

    return ET == 1;
}

bool PureMSA::ParityCheck(std::vector<char> &decoded_word){
    for (unsigned int idx=0; idx<BN.size(); idx++){
        double allR_L = BN[idx].received;
        for (unsigned int jdx=0; jdx<BN[idx].link.size(); jdx++){
            checkNode *linkedCN = &(CN[BN[idx].link[jdx]]);
            for (unsigned int cidx=0; cidx<linkedCN->link.size(); cidx++)
                if(linkedCN->link[cidx] == idx){
                    BN[idx].Q[jdx] = linkedCN->R[cidx];
                    allR_L += BN[idx].Q[jdx];
                    break;
                }
        }

        decoded_word[idx] = (allR_L > 0 ? 0:1);

    }
    matrixMultiplication(decoded_word,Hmatrix,parity_check);
    return oneCount(parity_check) == 0;
}