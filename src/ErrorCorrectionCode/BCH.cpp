#include "BCH.h"
#include "libMath.h"
#include "libDebug.h"
#include <cmath>
#include <cassert>
#include <cstdio>

BCH::BCH(std::map<std::string, std::string> &config) : ErrorCorrectionCodeBase(config){
	// default setting
	ECC_info = "BCH";
	BCH_order = 0;

	// deal config
	for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
		if      (pair_idx->first == "codeword_length")	 codeword_length = stoi(pair_idx->second);
		else if (pair_idx->first == "message_length")	 message_length = stoi(pair_idx->second);
		else WARRING(" key : %s is not find in BCH.cpp.",pair_idx->first.c_str());
	}

	assert(codeword_length != 0);
	if(matrix_src == "byAlgorithm"){
		if((codeword_length & (codeword_length+1)) != 0){
			ERROR("codeword_length %d is not valid.",codeword_length);
			exit(0);
		}
		while(codeword_length>>BCH_order != 0)
			BCH_order++;

		// build field
		std::vector<unsigned int> temp_v;
		creatFieldPoly(2,BCH_order,temp_v);
		field_poly = fieldPoly2Num(temp_v,2);

		DtoA.resize(codeword_length+1);
		AtoD.resize(codeword_length+1);
		DtoA[0] = 1;
		for(unsigned int idx=1; idx<DtoA.size();idx++){
			DtoA[idx] = DtoA[idx-1] << 1;
			if(DtoA[idx] > codeword_length)
				DtoA[idx] ^= field_poly;
			AtoD[DtoA[idx]] = idx%63;
			// INFO("%d %d",idx,DtoA[idx]);
		}

		// build valid_message_length
		std::vector<unsigned int> valid_message_length;
		std::vector<char> temp_char_v(codeword_length,0);
		unsigned int current_message_length = codeword_length;
		valid_message_length.push_back(current_message_length);
		for(unsigned int idx=1; current_message_length!=0 ;idx+=2){
			unsigned int target_idx = idx;
			for(unsigned int jdx=0; jdx<BCH_order; jdx++){
				temp_char_v[target_idx] = 1;
				target_idx = (target_idx*2) % (codeword_length);
			}
			current_message_length = codeword_length-oneCount(temp_char_v);
			if(valid_message_length[valid_message_length.size()-1]!=current_message_length)
				valid_message_length.push_back(current_message_length);
		}

		// check vaild message_length
		for(unsigned int idx=0; idx<valid_message_length.size()-1; idx++){
			if(valid_message_length[idx] == message_length)
				break;
			if(idx==valid_message_length.size()-2){
				printf("message_length %d is not valid.\n",message_length);
				printf("Here is all valid message_length, plz select one of these and try again.\n");
				printf("valid message_length :\n");
				for(unsigned int jdx=0; jdx<valid_message_length.size()-1; jdx++)
					printf("message_length : %d\n",valid_message_length[jdx]);
			}
		}

		// build G_poly
		temp_char_v.assign(codeword_length,0);
		current_message_length = codeword_length;
		unsigned int root_order_max = 0;
		for(unsigned int idx=1; current_message_length!=0 ;idx+=2){
			unsigned int target_idx = idx;
			root_order_max++;
			for(unsigned int jdx=0; jdx<BCH_order; jdx++){
				temp_char_v[target_idx] = 1;
				target_idx = (target_idx*2) % (codeword_length);
			}
			current_message_length = codeword_length-oneCount(temp_char_v);
			if(current_message_length==message_length)
				break;
		}
		G_poly.assign(1,1);
		std::vector<unsigned int> temp_G_poly(1,1);
		std::vector<unsigned int> alpha(2,1);
		for(unsigned int idx=0; idx<temp_char_v.size(); idx++){
			if(temp_char_v[idx]==0)
				continue;
			alpha[0] = DtoA[idx];
			fieldPolyMulti(temp_G_poly,alpha,G_poly,AtoD,DtoA);
			temp_G_poly.assign(G_poly.begin(), G_poly.end());
		}

		codeword_poly.resize(codeword_length);
		parity.resize(codeword_length - message_length);

		parity_check.resize(codeword_length - message_length);

		// build Gmatrix
		Gmatrix.resize(message_length);
		for(unsigned int idx=0; idx<message_length; idx++){
			Gmatrix[idx].assign(codeword_length,0);
			for(unsigned int jdx=0; jdx<G_poly.size(); jdx++)
				Gmatrix[idx][idx+jdx] = G_poly[jdx];
		}

		// build Hmatrix
		Hmatrix.resize(codeword_length);
		for(unsigned int idx=0; idx<codeword_length; idx++){
			Hmatrix[idx].assign(BCH_order*root_order_max,0);
			for(unsigned int root_idx=0; root_idx<root_order_max; root_idx++){
				unsigned int root = DtoA[(root_idx*2+1)*idx %(AtoD.size()-1)];
				for(unsigned int order_idx=0; order_idx<BCH_order; order_idx++)
					Hmatrix[idx][root_idx*BCH_order+BCH_order-1-order_idx] = (root >> order_idx) & 0x1;
			}
		}		

	}else{
		ERROR("BCH only support matrix_src=byAlgorithm.");
		exit(1);
	}

	// std::vector<std::vector<char> > check_matrix;
	// check_matrix.resize(Gmatrix.size());
	// for(unsigned int idx=0; idx<Gmatrix.size(); idx++)
	// 	check_matrix[idx].resize(Hmatrix[0].size());
	// printf("A\n");
	// matrixMultiplication(Gmatrix,Hmatrix,check_matrix);
	// printf("B\n");
	// for(unsigned int idx=0; idx<check_matrix.size(); idx++){
	// 	for(unsigned int jdx=0; jdx<check_matrix[0].size(); jdx++)
	// 		printf("%d",check_matrix[idx][jdx]);
	// 	printf("\n");
	// }
}

BCH::~BCH(){}

bool BCH::doEncode(std::vector<char> &message, std::vector<char> &codeword){

	codeword_poly.assign(codeword_length,0);
	for(unsigned int idx=0; idx<message.size(); idx++){
		codeword_poly[idx] = (unsigned int)message[idx];
	}
	fieldPolyMod(codeword_poly,G_poly,parity,AtoD,DtoA);
	for(unsigned int idx=0; idx<message.size(); idx++){
		codeword[idx] = message[idx];
	}
	for(unsigned int idx=0; idx<parity.size(); idx++){
		codeword[codeword_length-idx-1] = parity[parity.size()-idx-1];
	}
	return 0;
}

bool BCH::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
	ERROR("not finish");
	exit(1);
	return 0;
}
