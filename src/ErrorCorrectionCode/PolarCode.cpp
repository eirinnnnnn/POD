#include "PolarCode.h"
#include "libMath.h"
#include "libDebug.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>

PolarCode::PolarCode(std::map<std::string, std::string> config) : ErrorCorrectionCodeBase(config){
	// default setting
	ECC_info = "PolarCode";
	encoder_enable = false;
	decoder_enable = false;

	target_raw_BER = 0.5;
	list_size = 0;

	PolarCode_order = 0;
	
	// deal config
	for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
		if 		(pair_idx->first == "target_raw_BER")	 target_raw_BER = stod(pair_idx->second);
		else if (pair_idx->first == "list_size")	 	 list_size = stoi(pair_idx->second);
		else if (pair_idx->first == "codeword_length" && codeword_length == 0)			codeword_length = stoi(pair_idx->second);
		else if (pair_idx->first == "message_length" && message_length == 0)			message_length = stoi(pair_idx->second);
		else WARRING(" key : %s is not find in PolarCode.cpp.",pair_idx->first.c_str());
	}

	while((codeword_length>>PolarCode_order) != 0)
		PolarCode_order++;
	PolarCode_order--;
	INFO("PolarCode_order = %d",PolarCode_order);

	// assert 
	if((1<<PolarCode_order) != (int)codeword_length){
		ERROR("(1<<PolarCode_order) != codeword_length");	exit(1);	}
	if(target_raw_BER <= 0){
		ERROR("target_raw_BER < 0");	exit(1);	}
	if(list_size < 1){
		ERROR("list_size < 1");	exit(1);	}

	bhattacharyya.resize(codeword_length);
	bhattacharyya_order.assign(codeword_length,-1);
	relation_ship.resize(codeword_length);
	diverge_flag.resize(codeword_length);
	received_order.resize(codeword_length);
	Gmatrix.resize(message_length);
	if(matrix_src == "byAlgorithm"){
		setTargetRBER(target_raw_BER);
		buildRelationShip();
		encoder_enable = true;
		decoder_enable = true;
	}
	// SCL_memory[list][order][codeword_idx] resize
	list_enable.resize(list_size);
	path_metric.resize(list_size);
	extend_path_metric.resize(list_size*2);
	SCL_memory.resize(list_size);
	for(unsigned int list_idx=0; list_idx<SCL_memory.size(); list_idx++){
		SCL_memory[list_idx].resize(PolarCode_order+1);
		for(unsigned int order_idx=0; order_idx<SCL_memory[list_idx].size(); order_idx++)
			SCL_memory[list_idx][order_idx].resize(codeword_length);
	}
}

PolarCode::~PolarCode(){}

bool PolarCode::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
	if(!decoder_enable){
		ERROR("!decoder_enable");	exit(1);
	}
	// initial
	list_enable.assign(list_size,0);
	path_metric.assign(list_size,0);
	for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++)
		SCL_memory[0][PolarCode_order][received_order[codeword_idx]].node_value = received[codeword_idx]; // need to change to LLR domain
	list_enable[0] = 1;

	// INFO();
	// for(unsigned int idx=0; idx<SCL_memory[0][0].size(); idx++){
	// 	for(unsigned int jdx=0; jdx<SCL_memory[0].size(); jdx++)
	// 		printf("%lf ", SCL_memory[0][jdx][idx].node_value);
	// 	printf("\n");
	// }

	// compute SCL_memory & path_metric
	for(unsigned int decode_idx=0; decode_idx<received.size(); decode_idx++){
		// INFO("%d %d",decode_idx,diverge_flag[decode_idx]);
		for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++)
			if(list_enable[list_idx])
				computeNode(list_idx, 0, decode_idx);
		// printf("done ");
		if(diverge_flag[decode_idx])
			infoProcess(decode_idx);
		else
			frozenProcess(decode_idx);
		// printf("done ");
		if((decode_idx & 0x1) == 1)
			for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++)
				if(list_enable[list_idx])
					updateNode(list_idx, 0, decode_idx);
		// printf("done\n");
	}

	// check info, frozen, dynamic_forzen
	for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
		if(list_enable[list_idx])
			for(unsigned int decode_idx=0; decode_idx<received.size(); decode_idx++){
				int check_sum = 0;
				for(unsigned int check_idx=0; check_idx <relation_ship[decode_idx].size(); check_idx++)
					check_sum ^= SCL_memory[list_idx][0][relation_ship[decode_idx][check_idx]].node_HD;
				if(check_sum){
					list_enable[list_idx] = 0;
					break;
				}
			}
	}



	// for(unsigned int idx=0; idx<SCL_memory[0][0].size(); idx++){
	// 	for(unsigned int jdx=0; jdx<SCL_memory[0].size(); jdx++)
	// 		printf("%lf ", SCL_memory[0][jdx][idx].node_value);
	// 	printf("\n");
	// }

	// pick min path_metric
	double min_path_metric = 999999999;
	for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
		// printf("%2d %d %lf\n",list_idx, list_enable[list_idx],path_metric[list_idx]);
		if(list_enable[list_idx] && path_metric[list_idx] < min_path_metric){
			min_path_metric = path_metric[list_idx];
			for(unsigned int decode_idx=0; decode_idx<received.size(); decode_idx++){
				decoded_word[decode_idx] = SCL_memory[list_idx][PolarCode_order][received_order[decode_idx]].node_HD;
			}
		}
	}
	return 0;
}

bool PolarCode::setTargetRBER(double RBER){
	INFO("setTargetRBER %lf",RBER);
	// compute bhattacharyya based on log(-log) domain
	for(unsigned int idx=0; idx < codeword_length; idx++){
		bhattacharyya[idx] = log(-log(4*RBER*(1-RBER))/2);
		for(int order_idx=PolarCode_order-1; order_idx>=0; order_idx--){
			double tmp_double = bhattacharyya[idx];
			if((idx>>order_idx) & 0x1){
				bhattacharyya[idx] = log(2) + tmp_double;
			}else{
				if(tmp_double < -15){
					bhattacharyya[idx] = tmp_double*2;
				}else{
					tmp_double = exp(tmp_double)*-1;
					bhattacharyya[idx] = log(-tmp_double-log(2-exp(tmp_double)));
				}
			}

		}
	}
	// order bhattacharyya
	for(unsigned int idx=0; idx < codeword_length; idx++){
		double tmp_double = 999999999.9;
		double target_idx = -1;
		for(unsigned int sub_idx=0; sub_idx < codeword_length; sub_idx++){
			if(bhattacharyya_order[sub_idx] == -1 && bhattacharyya[sub_idx] < tmp_double){
				tmp_double = bhattacharyya[sub_idx];
				target_idx = sub_idx;
			}
		}
		bhattacharyya_order[target_idx] = idx;
	}
	return 0;
}

bool PolarCode::buildRelationShip(){
	INFO("buildRelationShip");
	// set relation_ship & Gmatrix & received_order
	for(unsigned int idx=0; idx < codeword_length; idx++){
		received_order[idx] = idx;
		INFO("index %d, order %d, bhattacharyya: %lf",idx, bhattacharyya_order[idx], bhattacharyya[idx]);
		if(bhattacharyya_order[idx] < int(codeword_length - message_length)){
			diverge_flag[idx] = 0;
			relation_ship[idx].assign(1,idx);
		}else{
			diverge_flag[idx] = 1;
			relation_ship[idx].resize(0);
		}
	}
	for(unsigned int idx=0,Gmatrix_idx=0; idx < codeword_length; idx++){
		if(relation_ship[idx].size() == 0){
			Gmatrix[Gmatrix_idx].resize(codeword_length);
			for(unsigned int sub_idx=0; sub_idx < codeword_length; sub_idx++){
				if((sub_idx | idx) == idx)
					Gmatrix[Gmatrix_idx][sub_idx] = 1;
				else
					Gmatrix[Gmatrix_idx][sub_idx] = 0;
			}
			Gmatrix_idx++;
		}
	}
	exit(0);
	return 0;
}

void PolarCode::computeNode(unsigned int list_idx, unsigned int order_idx, unsigned int decode_idx){
	if(order_idx == PolarCode_order)
		return;

	unsigned int pair_idx = decode_idx ^ (1<<order_idx);

	if(decode_idx < pair_idx){
		computeNode(list_idx,order_idx+1,pair_idx);
		computeNode(list_idx,order_idx+1,decode_idx);

		SCL_memory[list_idx][order_idx][decode_idx].node_value = SIGN(SCL_memory[list_idx][order_idx+1][decode_idx].node_value)
															   * SIGN(SCL_memory[list_idx][order_idx+1][pair_idx].node_value)
															   * MIN(ABS(SCL_memory[list_idx][order_idx+1][decode_idx].node_value),
															   	     ABS(SCL_memory[list_idx][order_idx+1][pair_idx].node_value));
	}else{
		SCL_memory[list_idx][order_idx][decode_idx].node_value = SCL_memory[list_idx][order_idx+1][decode_idx].node_value
															   + SCL_memory[list_idx][order_idx+1][pair_idx].node_value
															   * (SCL_memory[list_idx][order_idx][pair_idx].node_HD ? -1:1);
	}
}

void PolarCode::infoProcess(unsigned int decode_idx){
	extend_path_metric.assign(list_size*2,-1);
	int extend_list_count=0;

	// extend path
	for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
		if(list_enable[list_idx]){
			extend_list_count += 2;
			if(SCL_memory[list_idx][0][decode_idx].node_value < 0){
				extend_path_metric[list_idx*2  ] = path_metric[list_idx] + ABS(SCL_memory[list_idx][0][decode_idx].node_value);
				extend_path_metric[list_idx*2+1] = path_metric[list_idx];
			}else{
				extend_path_metric[list_idx*2  ] = path_metric[list_idx];
				extend_path_metric[list_idx*2+1] = path_metric[list_idx] + ABS(SCL_memory[list_idx][0][decode_idx].node_value);
			}
		}
	}

	// maintain path_metric.size() = list_size
	while(extend_list_count > list_size){
		double tmp_double = 0;
		int target_idx = -1;
		for(unsigned int idx=0; idx<extend_path_metric.size(); idx++){
			if(extend_path_metric[idx] > tmp_double){
				tmp_double = extend_path_metric[idx];
				target_idx = idx;
			}
		}
		extend_path_metric[target_idx] = -1;
		extend_list_count--;
	}
	// disable list_enable
	for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
		if(extend_path_metric[list_idx*2] == -1 && extend_path_metric[list_idx*2+1] == -1)
			list_enable[list_idx] = 0;
	}

	// move page
	for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
		if(extend_path_metric[list_idx*2] != -1 && extend_path_metric[list_idx*2+1] != -1){
			for(unsigned int target_idx=0; target_idx<list_enable.size(); target_idx++)
				if(list_enable[target_idx] == 0){
					clonePage(list_idx,target_idx);
					list_enable[target_idx] = 1;
					extend_path_metric[target_idx*2] = extend_path_metric[list_idx*2];
					extend_path_metric[list_idx*2] = -1;
					break;
				}
		}
	}

	// set path_metric
	for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
		if(list_enable[list_idx]){
			if(extend_path_metric[list_idx*2] == -1){
				path_metric[list_idx] = extend_path_metric[list_idx*2+1];
				SCL_memory[list_idx][0][decode_idx].node_HD = 1;
			}else{
				path_metric[list_idx] = extend_path_metric[list_idx*2];
				SCL_memory[list_idx][0][decode_idx].node_HD = 0;
			}
		}
	}	
}

void PolarCode::frozenProcess(unsigned int decode_idx){
	for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
		if(list_enable[list_idx]){
			int check_sum = 0;
			for(unsigned int relation_idx=0; relation_idx<relation_ship[decode_idx].size(); relation_idx++)
				if(relation_ship[decode_idx][relation_idx] != (int)decode_idx)
					check_sum ^= SCL_memory[list_idx][0][relation_ship[decode_idx][relation_idx]].node_HD;
			SCL_memory[list_idx][0][decode_idx].node_HD = check_sum;
			if(SIGN(SCL_memory[list_idx][0][decode_idx].node_value) * (0.5- SCL_memory[list_idx][0][decode_idx].node_HD) < 0)
				path_metric[list_idx] += ABS(SCL_memory[list_idx][0][decode_idx].node_value);
		}
	}
}

void PolarCode::updateNode(unsigned int list_idx, unsigned int order_idx, unsigned int decode_idx){
	if(order_idx == PolarCode_order)
		return;

	unsigned int pair_idx = decode_idx ^ (1<<order_idx);

	SCL_memory[list_idx][order_idx+1][decode_idx].node_HD = SCL_memory[list_idx][order_idx][decode_idx].node_HD;
	SCL_memory[list_idx][order_idx+1][pair_idx].node_HD   = SCL_memory[list_idx][order_idx][decode_idx].node_HD
														  ^ SCL_memory[list_idx][order_idx][pair_idx].node_HD;

	if((decode_idx & (1<<(order_idx+1))) != 0 ){
		updateNode(list_idx,order_idx+1,decode_idx);
		updateNode(list_idx,order_idx+1,pair_idx);
	}
}

void PolarCode::clonePage(unsigned int src_idx, unsigned int target_idx){
	for(unsigned int order_idx=0; order_idx<SCL_memory[src_idx].size(); order_idx++){
		for(unsigned int idx=0; idx<SCL_memory[src_idx][order_idx].size(); idx++){
			SCL_memory[target_idx][order_idx][idx].node_value = SCL_memory[src_idx][order_idx][idx].node_value;
			SCL_memory[target_idx][order_idx][idx].node_HD = SCL_memory[src_idx][order_idx][idx].node_HD;
		}
	}
}
