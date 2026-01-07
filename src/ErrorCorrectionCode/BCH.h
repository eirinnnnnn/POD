#ifndef _BCH_H_
#define _BCH_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

class BCH : public ErrorCorrectionCodeBase{
public:
    BCH(std::map<std::string, std::string> &config);
    virtual ~BCH();
    virtual bool doEncode(std::vector<char> &message, std::vector<char> &codeword);
    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

    std::vector<unsigned int> getAtoD() {return AtoD;   };
    std::vector<unsigned int> getDtoA() {return DtoA;   };
protected:
	// set by config

    // in BCH
    unsigned int BCH_order;
    unsigned int field_poly;
    std::vector<unsigned int> G_poly;
    std::vector<unsigned int> AtoD;
    std::vector<unsigned int> DtoA;


    std::vector<unsigned int> parity;
    std::vector<unsigned int> codeword_poly;


};


#endif