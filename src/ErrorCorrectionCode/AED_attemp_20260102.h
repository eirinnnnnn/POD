#ifndef _ADJUST_POLAR_DECODER_AED_20260102_H_
#define _ADJUST_POLAR_DECODER_AED_20260102_H_

#include <map>
#include <string>
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
    virtual double getCodeRate();
protected:
    // set by config
    std::string permutation_src;
    long permutation_random_seed;
    double target_raw_BER;
    std::string operationArray;
    unsigned int list_size;
    std::string bha_value_setting;
    bool OnlyInit;

    // === AED config ===
    bool use_AED;
    unsigned int aed_L;
    std::string automorphism_src;

    std::vector<std::vector<char> > polar_matrix;
    std::vector<double> bhattacharyya_value;
    unsigned int stage;

    std::vector<std::vector<char> > permutation_matrix;
    std::vector<unsigned int> received_order;

    // === AED precomputed orders ===
    std::vector<std::vector<unsigned int> > received_order_set;

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

    // === AED helpers ===
    bool decode_with_order_scl(const std::vector<double> &received,
                               const std::vector<unsigned int> &order,
                               std::vector<char> &decoded_word,
                               double &metric,
                               std::string &log);

    bool doDecode_AED(std::vector<double> &received,
                      std::vector<char> &decoded_word,
                      std::string &log);

    bool load_automorphism_set();
    bool parse_order_from_block(const std::vector<std::vector<unsigned int> > &raw,
                                unsigned int offset,
                                std::vector<unsigned int> &order);
};

#endif
