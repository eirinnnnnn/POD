#ifndef _POLAR_CODE_H_
#define _POLAR_CODE_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

// basic Successive Cancellation List (SCL) decode

typedef struct{
    double node_value;
    char node_HD; 
}SCL_node;

class PolarCode : public ErrorCorrectionCodeBase{
public:
    PolarCode(std::map<std::string, std::string> config);
    virtual ~PolarCode();
    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

    // setTargetRBER will build and order bhattacharyya 
    bool setTargetRBER(double RBER);
    // build relation_ship and made Gmatrix;
    virtual bool buildRelationShip();


protected:
	// set by config
	double target_raw_BER;
	int list_size;

	unsigned int PolarCode_order;
	std::vector<unsigned int> received_order;
	std::vector<std::vector<int> > relation_ship; // relation_ship[i].size() = 0 mean information bit.
												  // 						 = 1 mean frozen bit.
												  // 						 > 1 dynamic frozen bit. 
	std::vector<char> diverge_flag;	// 1: (information, frozen, dynamic frozen) can be deverge.
									// 0: (frozen, dynamic frozen) extend by expect bit.
	std::vector<double> bhattacharyya; // based on log(-log) domain
	std::vector<int> bhattacharyya_order; // 0: min

	std::vector<char> list_enable;
	std::vector<double> path_metric;
	std::vector<double> extend_path_metric;
	// extend_path_metric[list_idx*2] mean extend to 0; extend_path_metric[list_idx*2+1] mean extend to 1;
	std::vector<std::vector<std::vector<SCL_node> > > SCL_memory;
	// SCL_memory[list][order][codeword_idx], received side order = PolarCode_order
	

	void computeNode(unsigned int list_idx, unsigned int order_idx, unsigned int decode_idx);
	void infoProcess(unsigned int decode_idx);
	void frozenProcess(unsigned int decode_idx);
	void updateNode(unsigned int list_idx, unsigned int order_idx, unsigned int decode_idx);
	void clonePage(unsigned int src_idx, unsigned int target_idx);
};


#endif