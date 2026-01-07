#ifndef _POLAR_MULTI_KERNAL_H_
#define _POLAR_MULTI_KERNAL_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

typedef struct{
    double value;
    char HD;
    bool linked;
}node;

using page = std::vector<std::vector<node> >;

class AdjustPolarDecoder : public ErrorCorrectionCodeBase{
public:
    AdjustPolarDecoder(std::map<std::string, std::string> config);
    virtual ~AdjustPolarDecoder();

    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);
    // virtual bool doEncode(std::vector<char> &message, std::vector<char> &codeword);
    virtual double getCodeRate();
protected:
    // set by config
    std::string permutation_src;
    long permutation_random_seed;
    // std::string dynamic_frozen_process;
    double target_raw_BER;
    std::string operationArray;
    unsigned int list_size;
    std::string bha_value_setting;
    bool OnlyInit;

    std::vector<std::vector<char> > polar_matrix;
    std::vector<double> bhattacharyya_value;
    unsigned int stage;

    std::vector<std::vector<char> > permutation_matrix;
    std::vector<unsigned int> received_order;
    std::vector<char> diverge_flag;
    std::vector<std::vector<unsigned int> > relation_ship;
    std::vector<char> list_enable;
    std::vector<double> path_metric;
    std::vector<double> extend_path_metric;
    std::vector<page> SCL_mem;
    double basic_bhattacharyya_value;

    void do_node_value(unsigned int index, unsigned int level, page &page_memory);
    void update_node_HD(unsigned int index, unsigned int level, page &page_memory);

    void infoProcess(unsigned int decode_idx);
    void frozenProcess(unsigned int decode_idx);
    void clonePage(unsigned int src, unsigned int dst);
    
};






#endif
