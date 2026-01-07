#include "MLD.h"
#include "libMath.h"
#include <cmath>
#include <cassert>
#include <cstdio>

MLD::MLD(std::map<std::string, std::string> &config) : ErrorCorrectionCodeBase(config){
	// default setting
	// deal config

	// generator codeword space
	// assert(message_length <= 16);
	// if(message_length > 16)
	// 	exit(1);

	dumpHmatrix(".");
	dumpGmatrix(".");
	
	std::vector<char> tmp_message;
	tmp_message.assign(message_length,0);
	codeword_space.resize((int)pow(2.0,message_length));
	for(unsigned int i_idx=0; i_idx<codeword_space.size(); i_idx++){
		codeword_space[i_idx].resize(codeword_length);
		for(unsigned int message_idx=0; message_idx<message_length; message_idx++){
			tmp_message[message_idx] ^= 1;
			if(tmp_message[message_idx] == 1)
				break;
		}
		matrixMultiplication(tmp_message,Gmatrix,codeword_space[i_idx]);
	}

	from_one.resize(codeword_length);
	from_zero.resize(codeword_length);

}

MLD::~MLD(){}

bool MLD::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
	double tmp_distance;
	double min_distance = 9999999999.9;

	for(unsigned int codeword_idx=0; codeword_idx<codeword_length;codeword_idx++){
		from_one[codeword_idx]  = pow(received[codeword_idx]+1,2);
		from_zero[codeword_idx] = pow(received[codeword_idx]-1,2);
	}

	for(unsigned int space_idx=0; space_idx<codeword_space.size(); space_idx++){

		tmp_distance = 0;
		for(unsigned int codeword_idx=0; codeword_idx<codeword_length; codeword_idx++){
			tmp_distance += (codeword_space[space_idx][codeword_idx] ? from_one[codeword_idx]:from_zero[codeword_idx]);
		}
		if(tmp_distance < min_distance){
			min_distance = tmp_distance;
			decoded_word.assign(codeword_space[space_idx].begin(),codeword_space[space_idx].end());
		}
		if(min_distance == 0)
			break;
	}
	log = "0";


	return 0;
}