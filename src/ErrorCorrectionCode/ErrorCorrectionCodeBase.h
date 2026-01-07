#ifndef _ERROR_CORRECTION_CODE_BASE_
#define _ERROR_CORRECTION_CODE_BASE_
#include <map>
#include <string>
#include <vector>

class ErrorCorrectionCodeBase{
public:
    ErrorCorrectionCodeBase(std::map<std::string, std::string> &config);
    virtual ~ErrorCorrectionCodeBase();
    virtual bool doEncode(std::vector<char> &message, std::vector<char> &codeword);
    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);

    unsigned int getMessageLength()    {return message_length;     };
    unsigned int getCodewordLength()   {return codeword_length;    };
    virtual double getCodeRate()   {printf("message_length/codeword_length = %d/%d\n",message_length,codeword_length); return 1.0*message_length/codeword_length;    };

    void dumpGmatrix(std::string path);
    void dumpHmatrix(std::string path);

    std::vector<std::vector<char> > getGmatrix(){return Gmatrix;     };
    std::vector<std::vector<char> > getHmatrix(){return Hmatrix;     };

protected:
    std::string ECC_info;
    bool encoder_enable;
    bool decoder_enable;
    // set by config
    std::string matrix_src;
    std::string Gmatrix_path;
    std::string Hmatrix_path;
    bool dump_Gmatrix_enable;
    bool dump_Hmatrix_enable;

    std::vector<std::vector<char> > Gmatrix;
    std::vector<std::vector<char> > Hmatrix;
    unsigned int codeword_length;
    unsigned int transmitted_length;
    unsigned int message_length;

    std::vector<char> parity_check;
};


#endif
