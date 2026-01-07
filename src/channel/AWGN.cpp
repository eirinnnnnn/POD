#include "AWGN.h"
#include "libMath.h"
#include "libDebug.h"
#include <cmath>
#include <cstdio>
#include <cassert>
#include <string>

AWGN::AWGN(std::map<std::string, std::string> &config) : channelBase(config){
	seed = -1;
	var = 0;

	seed = std::stol(seed_string);
	assert(step_type == "SNR" || step_type == "var");
	if(step_type == "SNR"){
		SNR = start;
		var = (double)(pow(10,-SNR/10)/(code_rate)/2.0);
	}else if(step_type == "var"){
		var = start;
	}else{
		ERROR("step_type = %s\n",step_type.c_str());
		exit(0);
	}

}

AWGN::~AWGN(){}

bool AWGN::addNoise(std::vector<char> &codeword, std::vector<double> &received){
	assert(codeword.size() == received.size());
	// printf("SNR, var, code_rate, BER %e %e %e %e\n", SNR, var, code_rate, NormalDistributionCDF(0,1,var));
	for(unsigned int i_idx=0; i_idx<codeword.size(); i_idx++){
		double tmp_double = (codeword[i_idx] ? -1.0:1.0 )+ NormalDistribution(&seed)*sqrt(var);
		// double LLR = log(NormalDistributionPDF(tmp_double,1,var)/NormalDistributionPDF(tmp_double,-1,var));
		received[i_idx] = tmp_double;
	}
	return 0;
}

bool AWGN::nextChannel(){
	if((start+step-end)*(start-end) <= 0)
		return false;
	start += step;

	if(step_type == "SNR"){
		SNR = start;
		var = (double)(pow(10,-SNR/10)/(code_rate)/2.0);
	}else if(step_type == "var"){
		var = start;
	}
	return true;
}

bool AWGN::setSeedString(std::string in){
	seed = std::stol(in);
	return 0;
}

bool AWGN::setCodeRate(double in){
	printf("channel set code rate: %lf \n",in);
	code_rate = in;
	if(step_type == "SNR"){
		var = (double)(pow(10,-SNR/10)/(code_rate)/2.0);
	}
	return 0;
}

std::string AWGN::getSeedString(){
	char info[300];
	sprintf(info, "%ld",seed);
	std::string ans = info;
	return ans;
}

std::string AWGN::getChannelInfo(){
	char info[300];
	sprintf(info, "%s = %5.2lf, rawBER = %5.2lf",step_type.c_str(),start,getRawBER());
	std::string ans = info;
	return ans;
}

double AWGN::getRawBER(){
	return NormalDistributionCDF(0,1,var);
}