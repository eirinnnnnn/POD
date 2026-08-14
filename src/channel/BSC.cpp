#include "BSC.h"
#include "libMath.h"
#include "libDebug.h"
#include <cmath>
#include <cstdio>
#include <cassert>
#include <string>

BSC::BSC(std::map<std::string, std::string> &config) : channelBase(config){
	seed = -1;
	p = 0;
	step_mode = "additive";

	for(std::map<std::string, std::string>::iterator pairIdx=config.begin(); pairIdx!=config.end(); pairIdx++)
		if(pairIdx->first == "step_mode") step_mode = pairIdx->second;
	assert(step_mode == "additive" || step_mode == "multiplicative");

	seed = std::stol(seed_string);
	assert(step_type == "p");
	p = start;
	assert(p > 0 && p < 1);
}

BSC::~BSC(){}

bool BSC::addNoise(std::vector<char> &codeword, std::vector<double> &received){
	assert(codeword.size() == received.size());
	double llr_magnitude = log((1.0-p)/p);
	for(unsigned int i_idx=0; i_idx<codeword.size(); i_idx++){
		bool flip = ran0(&seed) < p;
		char received_bit = codeword[i_idx] ^ (flip ? 1:0);
		received[i_idx] = received_bit ? -llr_magnitude : llr_magnitude;
	}
	return 0;
}

bool BSC::nextChannel(){
	double next_start = (step_mode == "multiplicative") ? start*step : start+step;
	if((next_start-end)*(start-end) <= 0)
		return false;
	start = next_start;

	assert(step_type == "p");
	p = start;
	assert(p > 0 && p < 1);
	return true;
}

bool BSC::setSeedString(std::string in){
	seed = std::stol(in);
	return 0;
}

bool BSC::setCodeRate(double in){
	printf("channel set code rate: %lf \n",in);
	code_rate = in;
	return 0;
}

std::string BSC::getSeedString(){
	char info[300];
	sprintf(info, "%ld",seed);
	std::string ans = info;
	return ans;
}

std::string BSC::getChannelInfo(){
	char info[300];
	sprintf(info, "%s = %9.6lf, rawBER = %9.6lf",step_type.c_str(),start,getRawBER());
	std::string ans = info;
	return ans;
}

double BSC::getRawBER(){
	return p;
}
