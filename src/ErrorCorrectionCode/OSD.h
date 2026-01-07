#ifndef _OSD_H_
#define _OSD_H_
#include <map>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

class OSD : public ErrorCorrectionCodeBase{
public:
    OSD(std::map<std::string, std::string> &config);
    virtual ~OSD();
    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

protected:
	// set by config
	unsigned int OSD_order;

	std::vector<int> received_order; // -1: unsorted; 0: max abs(received); 1: second abs(received)
	std::vector<double> sorted_received;
	std::vector<std::vector<char> > OSD_Gmatrix;
	std::vector<char> tmp_message;
	std::vector<char> tmp_codeword;
	std::vector<unsigned int> flip_idx;

};


#endif