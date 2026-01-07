#ifndef _POLAR_MULTI_KERNAL_H_
#define _POLAR_MULTI_KERNAL_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

#include "KernalManager.h"

class PolarMultiKernal : public ErrorCorrectionCodeBase{
public:
    PolarMultiKernal(std::map<std::string, std::string> config);
    virtual ~PolarMultiKernal();

    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);
    // virtual bool doEncode(std::vector<char> &message, std::vector<char> &codeword);
    
    // setTargetRBER will build and order bhattacharyya 
    // bool setTargetRBER(double RBER);

    // build relation_ship and made Gmatrix;
    // virtual bool buildRelationShip();


protected:
    // set by config
    std::string permutation_src;
    std::string dynamic_frozen_process; 
    long permutation_random_seed;
    std::string kernal_string;
    std::string bha_value_setting;
    double target_raw_BER;
    unsigned int list_size;

    KernalManager *KM;
    std::vector<std::vector<char> > polar_matrix;
    unsigned int stage;

    std::vector<std::vector<char> > permutation_matrix;
    std::vector<unsigned int> received_order;
    std::vector<char> diverge_flag;
    std::vector<std::vector<unsigned int> > relation_ship;
    // relation_ship[info bit idx].size = 0
    // relation_ship[(dynamic) forzen bit idx] last value is (dynamic) forzen bit idx

    // mem
    std::vector<char> list_enable;
    std::vector<double> path_metric;
    std::vector<double> extend_path_metric;
    std::vector<page> SCL_mem;
    std::vector<char> extend_message;

    void infoProcess(unsigned int decode_idx);
    void frozenProcess(unsigned int decode_idx);
    void clonePage(unsigned int src, unsigned int dst);
    
};






#endif