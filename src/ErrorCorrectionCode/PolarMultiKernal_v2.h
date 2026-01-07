#ifndef _POLAR_MULTI_KERNAL_V2_H_
#define _POLAR_MULTI_KERNAL_V2_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

#include "KernalManager.h"

class PolarMultiKernal_v2 : public ErrorCorrectionCodeBase{
public:
    PolarMultiKernal_v2(std::map<std::string, std::string> config);
    virtual ~PolarMultiKernal_v2();

    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

protected:
    // set by config
    std::string permutation_src;
    std::string dynamic_frozen_process;
    std::string kernal_string;
    double target_raw_BER;
    unsigned int list_size;

    KernalManager *KM;
    std::vector<std::vector<char> > polar_matrix;
    unsigned int stage;




    std::vector<unsigned int> permutation_vertor;
    std::vector<unsigned int> Pp;
    std::vector<unsigned int> receive_order;
    std::vector<char> receive_valid;
    std::vector<char> diverge_flag;
    std::vector<std::vector<unsigned int> > relation_ship;

    std::vector<std::vector<char> > PolarGP;
    std::vector<std::vector<char> > PolarGPH;
    std::vector<std::vector<char> > PolarGPH_T;
    std::vector<unsigned int> tempV;
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