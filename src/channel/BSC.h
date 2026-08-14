#ifndef _BSC_H_
#define _BSC_H_
#include "channelBase.h"
#include <map>
#include <string>
#include <vector>

class BSC : public channelBase{
public:
	BSC(std::map<std::string, std::string> &config);
    ~BSC();
    virtual bool addNoise(std::vector<char> &codeword, std::vector<double> &received);
    virtual bool nextChannel();

	virtual bool setSeedString(std::string in);
	virtual bool setCodeRate(double in);
    virtual std::string getSeedString();
    virtual std::string getChannelInfo();
    virtual double getRawBER();

protected:
	long seed;
	double p;			// crossover probability
	std::string step_mode;	// "additive": start += step (default); "multiplicative": start *= step
};

#endif
