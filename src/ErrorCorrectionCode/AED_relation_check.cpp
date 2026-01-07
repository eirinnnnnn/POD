#include "AED_relation_check.h"
#include <limits>

bool AdjustPolarDecoderRelation::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
    // Fall back to base behavior if AED is disabled.
    if(!use_AED)
        return AdjustPolarDecoder::doDecode(received, decoded_word, log);

    if(received_order_set.empty()){
        use_AED = false;
        return AdjustPolarDecoder::doDecode(received, decoded_word, log);
    }

    unsigned int perm_limit = received_order_set.size();
    bool found_any = false;
    double best_metric = std::numeric_limits<double>::max();
    std::vector<char> best_word(codeword_length, 0);

    // For each permutation, run full SCL (relation check is inside decode_with_order_scl).
    for(unsigned int perm_idx=0; perm_idx<perm_limit; perm_idx++){
        std::vector<char> candidate(codeword_length, 0);
        double metric = 0.0;
        if(!decode_with_order_scl(received, received_order_set[perm_idx], candidate, metric, log))
            continue;

        if(metric < best_metric){
            best_metric = metric;
            best_word.swap(candidate);
            found_any = true;
        }
    }

    if(found_any){
        decoded_word.swap(best_word);
        return true;
    }

    // If nothing survived, fall back to base SCL on baseline permutation.
    use_AED = false;
    return AdjustPolarDecoder::doDecode(received, decoded_word, log);
}
