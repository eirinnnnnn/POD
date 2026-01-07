#include "channelBase.h"
#include <cassert>
channelBase::channelBase(std::map<std::string, std::string> &config){
	step_type = "SNR";
	start = 0;
	step = 0;
	end = 0;
	seed_string = "";
	code_rate = 1;
	
	raw_BER = 0;
	SNR = 0;

	for(std::map<std::string, std::string>::iterator pairIdx=config.begin(); pairIdx!=config.end(); pairIdx++) {  
			if	 	(pairIdx->first == "step_type" )	step_type	= pairIdx->second;
			else if (pairIdx->first == "start")			start 		= std::stod(pairIdx->second);
			else if (pairIdx->first == "step")			step 		= std::stod(pairIdx->second);
			else if (pairIdx->first == "end")			end 		= std::stod(pairIdx->second);
			else if (pairIdx->first == "seed_string")	seed_string = pairIdx->second;
			else if (pairIdx->first == "code_rate")		code_rate	= std::stod(pairIdx->second);
	}
}

channelBase::~channelBase(){}

bool channelBase::addNoise(std::vector<char> &codeword, std::vector<double> &received){
	assert(codeword.size() == received.size());
	for(unsigned int i_idx=0; i_idx<codeword.size(); i_idx++)
		received[i_idx] = codeword[i_idx] ? -1.0:1.0;
	return 0;
};
