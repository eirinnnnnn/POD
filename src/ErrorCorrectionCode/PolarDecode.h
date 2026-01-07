#ifndef _POLAR_DECODE_H_
#define _POLAR_DECODE_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"
#include "PolarCode.h"

// basic Successive Cancellation List (SCL) decode

class PolarDecode : public PolarCode{
public:
    PolarDecode(std::map<std::string, std::string> config);
    virtual ~PolarDecode();
    // virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

    // build relation_ship and made Gmatrix;
    virtual bool buildRelationShip();
protected:
    // set by config
    std::string permutation_src;
    std::string dynamic_frozen_process;
    std::vector<std::vector<char> > permutation_matrix;
    std::vector<unsigned int> permutation_array;
    long permutation_random_seed;

    // std::vector<double> internal_received;
    // std::vector<char> internal_decoded_word;

};


#endif