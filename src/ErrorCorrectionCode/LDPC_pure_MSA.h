#ifndef _LDPC_PURE_MSA_H_
#define _LDPC_PURE_MSA_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

typedef struct{
    double received;
    std::vector<unsigned int> link;
    std::vector<double> Q;
}bitNode;

typedef struct{
    char allSign;
    std::vector<unsigned int> link;
    std::vector<double> R;
}checkNode;




class PureMSA : public ErrorCorrectionCodeBase{
public:
    PureMSA(std::map<std::string, std::string> &config);
    virtual ~PureMSA();
    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

protected:
    std::vector<bitNode> BN;
    std::vector<checkNode> CN;
    void BNP();
    char CNP_MSA();
    char CNP_SPA();
    bool ParityCheck(std::vector<char> &decoded_word);

    void printf_BN();
    void printf_CN();

    unsigned Maxiter;
};

#endif