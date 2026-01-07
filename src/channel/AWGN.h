#ifndef _AWGN_H_
#define _AWGN_H_
#include "channelBase.h"
#include <map>
#include <string>
#include <vector>

class AWGN : public channelBase{
public:
	AWGN(std::map<std::string, std::string> &config);
    ~AWGN();
    virtual bool addNoise(std::vector<char> &codeword, std::vector<double> &received);
    virtual bool nextChannel();
	
	virtual bool setSeedString(std::string in);
	virtual bool setCodeRate(double in);
    virtual std::string getSeedString();
    virtual std::string getChannelInfo();
    virtual double getRawBER();
    

protected:
	long seed;
	double var;
};

#endif