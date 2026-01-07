#include "OSD.h"
#include "libMath.h"
#include "libDebug.h"
#include <cmath>
#include <cassert>
#include <cstdio>

OSD::OSD(std::map<std::string, std::string> &config) : ErrorCorrectionCodeBase(config){
	// default setting
	OSD_order = 0;

	// deal config
	for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
		if      (pair_idx->first == "OSD_order")	 OSD_order = stoi(pair_idx->second);
		else WARRING(" key : %s is not find in OSD.cpp.",pair_idx->first.c_str());
	}

	received_order.resize(codeword_length);
	sorted_received.resize(codeword_length);
	tmp_message.resize(message_length);
	tmp_codeword.resize(codeword_length);
	flip_idx.resize(OSD_order);
	OSD_Gmatrix.resize(message_length);
	for(unsigned int idx=0; idx<message_length;idx++)
		OSD_Gmatrix[idx].resize(codeword_length);
}

OSD::~OSD(){}

bool OSD::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){

	// initial
	received_order.assign(codeword_length,-1);

	// received, OSD_Gmatrix sorting
	for(unsigned int loop=0; loop<codeword_length; loop++){
		double max_received = 0;
		unsigned int target_idx = -1;
		for(unsigned int idx=0; idx<codeword_length; idx++)
			if(received_order[idx] == -1 && max_received < ABS(received[idx])){
				target_idx = idx;
				max_received = ABS(received[idx]);
			}
		received_order[target_idx] = loop;
		for(unsigned int idx=0; idx<message_length; idx++)
			OSD_Gmatrix[idx][loop] = Gmatrix[idx][target_idx];
	}

	GaussianJordanElimination(OSD_Gmatrix,0);

	// make OSD_Gmatrix be systematic
	for(unsigned int row_idx=0; row_idx<OSD_Gmatrix.size(); row_idx++){
		if(OSD_Gmatrix[row_idx][row_idx])
			continue;
		unsigned int target_idx = -1;
		for(unsigned int col_idx = row_idx+1; col_idx<OSD_Gmatrix[0].size(); col_idx++)
			if(OSD_Gmatrix[row_idx][col_idx]){
				target_idx = col_idx;
				break;
			}
		// switch target_idx to row_idx
		for(unsigned int idx=0; idx<codeword_length; idx++){
			if(received_order[idx] == (int)target_idx)
				received_order[idx] = row_idx;
			else if(received_order[idx] == (int)row_idx)
				received_order[idx] = target_idx;
		}
		// unsigned int tmpInt = received_order[row_idx];
		// received_order[row_idx] = received_order[target_idx];
		// received_order[target_idx] = tmpInt;

		for(unsigned int sub_row_idx=0; sub_row_idx<OSD_Gmatrix.size(); sub_row_idx++){
			OSD_Gmatrix[sub_row_idx][target_idx] = OSD_Gmatrix[sub_row_idx][row_idx];
			OSD_Gmatrix[sub_row_idx][row_idx] = (sub_row_idx == row_idx ? 1:0);
		}
	}

	// get sorted_received
	for(unsigned int idx=0; idx<codeword_length; idx++)
		sorted_received[received_order[idx]] = received[idx];

	// OSD_order = 0; basic computing
	for(unsigned int idx=0; idx<message_length; idx++)
		tmp_message[idx] = (sorted_received[idx]<0 ? 1:0);
	matrixMultiplication(tmp_message,OSD_Gmatrix,tmp_codeword);
	double final_distance = 0;
	for(unsigned int idx=0; idx<codeword_length; idx++)
		final_distance += pow((tmp_codeword[idx] ? -1:1)-sorted_received[idx],2);
	for(unsigned int idx=0; idx<codeword_length; idx++)
		decoded_word[idx] = tmp_codeword[received_order[idx]];
	
	for(unsigned int order_idx=0; order_idx<OSD_order; order_idx++){
		// initial
		for(unsigned int idx=0; idx<=order_idx; idx++)
			flip_idx[idx] = idx;
		while(flip_idx[order_idx] != message_length){
			for(unsigned int idx=0; idx<=order_idx; idx++)
				tmp_message[flip_idx[idx]] ^= 1;

			matrixMultiplication(tmp_message,OSD_Gmatrix,tmp_codeword);
			double tmp_distance = 0;
			for(unsigned int idx=0; idx<codeword_length; idx++)
				tmp_distance += pow((tmp_codeword[idx] ? -1:1)-sorted_received[idx],2);
			if(tmp_distance < final_distance){
				final_distance = tmp_distance;
				for(unsigned int idx=0; idx<codeword_length; idx++)
					decoded_word[idx] = tmp_codeword[received_order[idx]];
			}

			for(unsigned int idx=0; idx<=order_idx; idx++)
				tmp_message[flip_idx[idx]] ^= 1;

			flip_idx[0]++;
			for(unsigned int idx=0; idx<order_idx; idx++)
				if(flip_idx[idx] == flip_idx[idx+1]){
					flip_idx[idx] = idx;
					flip_idx[idx+1]++;
				}
		}
	}
	return 0;
}
