#ifndef _MLD_H_
#define _MLD_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

class MLD : public ErrorCorrectionCodeBase{
public:
    MLD(std::map<std::string, std::string> &config);
    virtual ~MLD();
    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

protected:
    std::vector<std::vector<char> > codeword_space;
    std::vector<double> from_one;
    std::vector<double> from_zero;
};

#endif