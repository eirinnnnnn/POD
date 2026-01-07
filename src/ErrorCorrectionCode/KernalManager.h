#ifndef _KERNAL_MANAGER_H_
#define _KERNAL_MANAGER_H_
#include <string>
#include <vector>

typedef struct{
    double value;
    char HD; 
}node;

using page = std::vector<std::vector<node> >;

class KernalManager{
public:
    KernalManager(std::string kernal_type_setting);
    virtual ~KernalManager();

    void set_target_BER(double BER);
    void build_Gmatrix(std::vector<std::vector<char> > &matrix);

    unsigned int get_index(std::vector<unsigned int> index);
    std::vector<unsigned int> get_index_array(unsigned int index);
    unsigned int get_stage(){
        return stage;                       };
    unsigned int get_codeword_length(){
        return codeword_length;             };

    void   set_basic_bhattacharyya_value(double sbhv){
        basic_bhattacharyya_value = sbhv;   };
    double get_basic_bhattacharyya_value(){
        return basic_bhattacharyya_value;   };

    void do_node_value(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    void do_node_value(unsigned int index, unsigned int level, page &page_memory){
        do_node_value(get_index_array(index),level,page_memory);    };

    void update_node_HD(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    void update_node_HD(unsigned int index, unsigned int level, page &page_memory){
        update_node_HD(get_index_array(index),level,page_memory);   };    

    double get_bhattacharyya_value(std::vector<unsigned int> index, unsigned int level);
    double get_bhattacharyya_value(unsigned int index, unsigned int level){
        return get_bhattacharyya_value(get_index_array(index),level);   };
    
    std::vector<double> bhattacharyya_value;
    std::vector<unsigned int> bhattacharyya_order;
    std::vector<std::string> bha_value_type;
    
protected:
    std::vector<std::string> kernal_type_array;
    unsigned int codeword_length;
    unsigned int stage;
    std::vector<unsigned int> weight;

    double target_BER;
    double basic_bhattacharyya_value;
    double erase_bhattacharyya_value;
    double know_bhattacharyya_value;
    
    

    // kernal unit
    void dnv_PK23(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    void unh_PK23(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    double bhv_PK23(std::vector<unsigned int> index, unsigned int level);

    void dnv_PK753(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    void unh_PK753(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    double bhv_PK753(std::vector<unsigned int> index, unsigned int level);

    void dnv_PK427(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    void unh_PK427(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    double bhv_PK427(std::vector<unsigned int> index, unsigned int level);

    void dnv_PK657(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    void unh_PK657(std::vector<unsigned int> index, unsigned int level, page &page_memory);
    double bhv_PK657(std::vector<unsigned int> index, unsigned int level);

};



#endif