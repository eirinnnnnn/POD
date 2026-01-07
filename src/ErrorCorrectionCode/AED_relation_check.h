#ifndef _ADJUST_POLAR_DECODER_RELATION_CHECK_H_
#define _ADJUST_POLAR_DECODER_RELATION_CHECK_H_

#include "AED_attemp_20260102.h"

// Variant decoder: when AED is enabled, selection across permutations
// is done only with the internal relation check (no external H^T).
class AdjustPolarDecoderRelation : public AdjustPolarDecoder{
public:
    AdjustPolarDecoderRelation(std::map<std::string, std::string> config)
    : AdjustPolarDecoder(config){}

    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);
};

#endif
