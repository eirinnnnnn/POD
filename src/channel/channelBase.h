#ifndef _CHANNEL_BASE_H_
#define _CHANNEL_BASE_H_
#include <map>
#include <string>
#include <vector>

class channelBase{
public:
	channelBase(std::map<std::string, std::string> &config);
    ~channelBase();
    virtual bool addNoise(std::vector<char> &codeword, std::vector<double> &received);
    virtual bool nextChannel()						{return false;				};
	
	virtual bool setSeedString(std::string in)	{seed_string = in;		return 0;};
	virtual bool setCodeRate(double in)			{code_rate = in;		return 0;};
    virtual std::string getSeedString()			{return seed_string;	};
    virtual std::string getChannelInfo()		{return "";	};
    virtual double getRawBER()				{return raw_BER;		};
    // virtual double get_raw_SNR()					{return SNR;			};

protected:
	std::string step_type;
	double start,step,end;
	std::string seed_string;
	double code_rate;

	double raw_BER;
	double SNR;
};

#endif