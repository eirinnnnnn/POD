#include "KernalManager.h"
#include "libDebug.h"
#include "libMath.h"
#include <cmath>
// #include <cassert>
// #include <cstdio>
// #include <vector>


KernalManager::KernalManager(std::string kernal_type_setting){
    kernal_type_array.assign(0,"");
    codeword_length = 1;
    stage = 0;
    weight.assign(0,1);

    while(true){
        // analysis kernal_type_setting
        int comma_index = kernal_type_setting.find(",");
        if(comma_index == -1){
            codeword_length *= kernal_type_setting.size();
            kernal_type_array.push_back(kernal_type_setting);
            break;
        }else{
            codeword_length *= comma_index;
            kernal_type_array.push_back(kernal_type_setting.substr(0,comma_index));
        }
        kernal_type_setting = kernal_type_setting.substr(comma_index+1,kernal_type_setting.size());
    }

    bhattacharyya_value.assign(codeword_length,0);
    bhattacharyya_order.assign(codeword_length,0);
    bha_value_type.assign(codeword_length,"?");

    stage = kernal_type_array.size();
    weight.assign(stage,1);

    for(unsigned int idx=0; idx<stage; idx++){
        for(unsigned int sub_idx=idx+1; sub_idx<stage; sub_idx++)
            weight[sub_idx] *= kernal_type_array[idx].size();
    }


    // in log(-log(Z)) domain, log(-log(4*target_BER*(1-target_BER))/2);
    erase_bhattacharyya_value = log(-log(4*0.5*(1-0.5))/2);   
    know_bhattacharyya_value  = 99999999999999999;  
    set_target_BER(0.1);

    // debug info
    if(DEBUG_MODE){
        printf("--KernalManager--\n");
        printf("codeword_length %d, stage %d\n",codeword_length, stage);
        for(unsigned int idx=0; idx<stage; idx++){
            printf("kernal %s, weight %d\n",kernal_type_array[idx].c_str(),weight[idx]);
        }
        printf("--KernalManager--\n");
    }

}

KernalManager::~KernalManager(){}

void KernalManager::set_target_BER(double BER){

    target_BER = BER;
    basic_bhattacharyya_value = log(-log(4*target_BER*(1-target_BER))/2);
    for(unsigned int idx=0; idx<codeword_length; idx++){
        bhattacharyya_value[idx] = get_bhattacharyya_value(idx,0);
    }
    
    // sorting, max = 0
    bhattacharyya_order.assign(codeword_length,0);
    double tmp_double;
    unsigned int site;
    for(unsigned int idx=0; idx<codeword_length; idx++){
        tmp_double = 9999999999;
        for(unsigned int jdx=0; jdx<codeword_length; jdx++){
            if(bhattacharyya_order[jdx] == 0 && bhattacharyya_value[jdx]<tmp_double){
                site = jdx;
                tmp_double = bhattacharyya_value[jdx];
            }
        }
        bhattacharyya_order[site] = codeword_length-idx-1;
    }
}

void KernalManager::build_Gmatrix(std::vector<std::vector<char> > &matrix){
    char g23[2][2] = {{1,0},
                      {1,1}};
    char g753[3][3] = {{1,1,1},
                       {1,0,1},
                       {0,1,1}};
    char g427[3][3] = {{1,0,0},
                       {0,1,0},
                       {1,1,1}};
    char g657[3][3] = {{1,1,0},
                       {1,0,1},
                       {1,1,1}};

    matrix.resize(codeword_length);
    for(unsigned int idx=0; idx<codeword_length; idx++){
        std::vector<unsigned int> idx_array = get_index_array(idx);
        matrix[idx].resize(codeword_length);
        for(unsigned int jdx=0; jdx<codeword_length; jdx++){
            std::vector<unsigned int> jdx_array = get_index_array(jdx);
            matrix[idx][jdx] = 1;
            for(unsigned int satge_idx=0; satge_idx<stage; satge_idx++){
                if      (kernal_type_array[satge_idx] == "23"){
                    matrix[idx][jdx] &=  g23[idx_array[satge_idx]][jdx_array[satge_idx]];
                }else if(kernal_type_array[satge_idx] == "753"){
                    matrix[idx][jdx] &= g753[idx_array[satge_idx]][jdx_array[satge_idx]];
                }else if(kernal_type_array[satge_idx] == "427"){
                    matrix[idx][jdx] &= g427[idx_array[satge_idx]][jdx_array[satge_idx]];
                }else if(kernal_type_array[satge_idx] == "657"){
                    matrix[idx][jdx] &= g657[idx_array[satge_idx]][jdx_array[satge_idx]];
                }else
                    ERROR("build_Gmatrix fail. %s",kernal_type_array[satge_idx].c_str());
            }
        }
    }
}

std::vector<unsigned int> KernalManager::get_index_array(unsigned int index){
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(stage,0);
    for(unsigned int idx=0; idx<stage; idx++){
        tmp_index[idx] = index % kernal_type_array[idx].size();
        index = index / kernal_type_array[idx].size();
    }
    return tmp_index; 
}

unsigned int KernalManager::get_index(std::vector<unsigned int> index){
    unsigned int ans=0;
    for(unsigned int idx=0; idx<stage; idx++)
        ans += (index[idx]*weight[idx]);
    return ans;
}

void KernalManager::do_node_value(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    INFO("dnv idx, level =  %3d, %d", get_index(index), level);
    // level end
    if(level == index.size())
        return; // return LLR

    // call kernal_unit
    if      (kernal_type_array[level] == "23"){
        dnv_PK23(index,level,page_memory);
    }else if(kernal_type_array[level] == "753"){
        dnv_PK753(index,level,page_memory);
    }else if(kernal_type_array[level] == "427"){
        dnv_PK427(index,level,page_memory);
    }else if(kernal_type_array[level] == "657"){
        dnv_PK657(index,level,page_memory);
    }else
        ERROR("kernal fail. %s",kernal_type_array[level].c_str());
}

void KernalManager::update_node_HD(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    INFO("unh idx, level =  %3d, %d", get_index(index), level);
    // level end
    if(level == index.size())
        return;

    // call kernal_unit
    if      (kernal_type_array[level] == "23"){
        unh_PK23(index,level,page_memory);
    }else if(kernal_type_array[level] == "753"){
        unh_PK753(index,level,page_memory);
    }else if(kernal_type_array[level] == "427"){
        unh_PK427(index,level,page_memory);
    }else if(kernal_type_array[level] == "657"){
        unh_PK657(index,level,page_memory);
    }else
        ERROR("kernal fail. %s",kernal_type_array[level].c_str());
}

double KernalManager::get_bhattacharyya_value(std::vector<unsigned int> index, unsigned int level){
    // level end
    if(level == index.size()){
        int ori_idx = index[level-1];
        for (int idx = index.size()-2; idx>=0; idx--){
            ori_idx = ori_idx * kernal_type_array[idx].size() + index[idx];
        }
        if (bha_value_type[ori_idx] == "0")
            return erase_bhattacharyya_value;
        else if (bha_value_type[ori_idx] == "1")
            return know_bhattacharyya_value;
        else if (bha_value_type[ori_idx] == "?")
            return basic_bhattacharyya_value;
        else
            ERROR("bha_value_type[%d]=%s", ori_idx, bha_value_type[ori_idx].c_str());
    }

    // call kernal_unit
    double return_value;
    if      (kernal_type_array[level] == "23"){
        return_value = bhv_PK23(index,level);
    }else if(kernal_type_array[level] == "753"){
        return_value = bhv_PK753(index,level);
    }else if(kernal_type_array[level] == "427"){
        return_value = bhv_PK427(index,level);
    }else if(kernal_type_array[level] == "657"){
        return_value = bhv_PK657(index,level);
    }else
        ERROR("kernal fail. %s",kernal_type_array[level].c_str());
    return return_value;
}

// -----------------------------
// kernal unit
// -----------------------------

void KernalManager::dnv_PK23(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    INFO("kernal PK23 level %d",level);
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1;
    switch(index[level]){
        case 0:
            // need to get source
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);

            page_memory[level][site_0].value = SIGN(page_memory[level+1][site_0].value) 
                                             * SIGN(page_memory[level+1][site_1].value)
                                             * MIN(ABS(page_memory[level+1][site_0].value),
                                                   ABS(page_memory[level+1][site_1].value));
            break;
        case 1:
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);

            page_memory[level][site_1].value =  page_memory[level+1][site_1].value
                                             + (page_memory[level][site_0].HD ? -1:1)
                                             *  page_memory[level+1][site_0].value;
            break;
        default:
            ERROR("kernal PK23 index with %d\n",index[level]);
            break;
    }
}

void KernalManager::unh_PK23(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    // not ready update
    if(index[level] != 1)
        return; 

    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1;

    tmp_index[level] = 0;   site_0 = get_index(tmp_index);
    tmp_index[level] = 1;   site_1 = get_index(tmp_index);

    page_memory[level+1][site_0].HD = page_memory[level][site_0].HD 
                                    ^ page_memory[level][site_1].HD;
    page_memory[level+1][site_1].HD = page_memory[level][site_1].HD;

    tmp_index[level] = 0;   update_node_HD(tmp_index, level+1, page_memory);
    tmp_index[level] = 1;   update_node_HD(tmp_index, level+1, page_memory);
}

double KernalManager::bhv_PK23(std::vector<unsigned int> index, unsigned int level){
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    double Z0, Z1;

    tmp_index[level] = 0;   Z0 = get_bhattacharyya_value(tmp_index, level+1);
    tmp_index[level] = 1;   Z1 = get_bhattacharyya_value(tmp_index, level+1);
    
    Z0 = exp(-exp(Z0));
    Z1 = exp(-exp(Z1));

    switch(index[level]){
        case 0:
            // Z0+Z1-Z0Z1
            return log(-log(Z0+Z1-Z0*Z1));
            break;
        case 1:
            // Z0Z1
            return log(-log(Z0*Z1));
            break;
        default:
            ERROR("bhv_PK23 index with %d\n",index[level]);
            break;
    }   
    return 0;
}

void KernalManager::dnv_PK753(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    INFO("kernal PK753 level %d",level);
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1, site_2;
    switch(index[level]){
        case 0:
            // need to get source
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);

            // u0 = x0x1x2
            page_memory[level][site_0].value = SIGN(page_memory[level+1][site_0].value) 
                                             * SIGN(page_memory[level+1][site_1].value)
                                             * SIGN(page_memory[level+1][site_2].value)
                                             * MIN(ABS(page_memory[level+1][site_0].value),
                                                   MIN(ABS(page_memory[level+1][site_1].value),
                                                       ABS(page_memory[level+1][site_2].value)));

            break;
        case 1:
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);

            // u1 = x1x2 or u0x0
            page_memory[level][site_1].value = SIGN(page_memory[level+1][site_1].value) 
                                             * SIGN(page_memory[level+1][site_2].value)
                                             * MIN(ABS(page_memory[level+1][site_1].value),
                                                   ABS(page_memory[level+1][site_2].value))
                                             + (page_memory[level][site_0].HD ? -1:1)
                                             *  page_memory[level+1][site_0].value;
            break;
        case 2:
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);

            // u2 = u0x1 or u0u1x2
            page_memory[level][site_2].value = (page_memory[level][site_0].HD ? -1:1)
                                             *  page_memory[level+1][site_1].value
                                             + (page_memory[level][site_0].HD ^ page_memory[level][site_1].HD ? -1:1)
                                             *  page_memory[level+1][site_2].value;

            break;
        default:
            ERROR("kernal PK753 index with %d\n",index[level]);
            break;
    }
}

void KernalManager::unh_PK753(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    // not ready update
    if(index[level] != 2)
        return; 

    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1, site_2;

    tmp_index[level] = 0;   site_0 = get_index(tmp_index);
    tmp_index[level] = 1;   site_1 = get_index(tmp_index);
    tmp_index[level] = 2;   site_2 = get_index(tmp_index);

    page_memory[level+1][site_0].HD = page_memory[level][site_0].HD 
                                    ^ page_memory[level][site_1].HD;
    page_memory[level+1][site_1].HD = page_memory[level][site_0].HD 
                                    ^ page_memory[level][site_2].HD;
    page_memory[level+1][site_2].HD = page_memory[level][site_0].HD
                                    ^ page_memory[level][site_1].HD
                                    ^ page_memory[level][site_2].HD;

    tmp_index[level] = 0;   update_node_HD(tmp_index, level+1, page_memory);
    tmp_index[level] = 1;   update_node_HD(tmp_index, level+1, page_memory);
    tmp_index[level] = 2;   update_node_HD(tmp_index, level+1, page_memory);
}

double KernalManager::bhv_PK753(std::vector<unsigned int> index, unsigned int level){
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    double Z0, Z1, Z2;

    tmp_index[level] = 0;   Z0 = get_bhattacharyya_value(tmp_index, level+1);
    tmp_index[level] = 1;   Z1 = get_bhattacharyya_value(tmp_index, level+1);
    tmp_index[level] = 2;   Z2 = get_bhattacharyya_value(tmp_index, level+1);
    
    Z0 = exp(-exp(Z0));
    Z1 = exp(-exp(Z1));
    Z2 = exp(-exp(Z2));

    switch(index[level]){
        case 0:
            // Z0+Z1+Z2-Z0Z1-Z1Z2-Z2Z0+Z0Z1Z2
            return log(-log(Z0+Z1+Z2-Z0*Z1-Z1*Z2-Z2*Z0+Z0*Z1*Z2));
            break;
        case 1:
            // Z0Z1 + Z0Z2 - Z0Z1Z2
            return log(-log(Z0*Z1 + Z0*Z2 - Z0*Z1*Z2));
            break;
        case 2:
            // Z1*Z2
            return log(-log(Z1*Z2));
            break;
        default:
            ERROR("bhv_PK753 index with %d\n",index[level]);
            break;
    }   
    return 0;
}

void KernalManager::dnv_PK427(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1, site_2;
    double sub_node = 0;

    switch(index[level]){
        case 0:
            // need to get source
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);

            // u0 = x0x2
            page_memory[level][site_0].value = SIGN(page_memory[level+1][site_0].value) 
                                             * SIGN(page_memory[level+1][site_2].value)
                                             * MIN(ABS(page_memory[level+1][site_0].value),
                                                   ABS(page_memory[level+1][site_2].value));
            break;
        case 1:
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);

            // u1 = x1 (x2 or u0x0)
            sub_node = page_memory[level+1][site_2].value
                    + (page_memory[level][site_0].HD ? -1:1) * page_memory[level+1][site_0].value;

            page_memory[level][site_1].value = SIGN(page_memory[level+1][site_1].value) 
                                             * SIGN(sub_node)
                                             * MIN(ABS(page_memory[level+1][site_1].value),
                                                   ABS(sub_node));

            break;
        case 2:
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);

            // u2 = u0x0 or u1x1 or  x2 
            page_memory[level][site_2].value = (page_memory[level][site_0].HD ? -1:1) * page_memory[level+1][site_0].value
                                             + (page_memory[level][site_1].HD ? -1:1) * page_memory[level+1][site_1].value
                                             +                                            page_memory[level+1][site_2].value;

            break;
        default:
            ERROR("dnv_PK427 index with %d\n",index[level]);
            break;
    }
}

void KernalManager::unh_PK427(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    // not ready update
    if(index[level] != 2)
        return; 

    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1, site_2;

    tmp_index[level] = 0;   site_0 = get_index(tmp_index);
    tmp_index[level] = 1;   site_1 = get_index(tmp_index);
    tmp_index[level] = 2;   site_2 = get_index(tmp_index);

    page_memory[level+1][site_0].HD = page_memory[level][site_0].HD 
                                    ^ page_memory[level][site_2].HD;
    page_memory[level+1][site_1].HD = page_memory[level][site_1].HD 
                                    ^ page_memory[level][site_2].HD;
    page_memory[level+1][site_2].HD = page_memory[level][site_2].HD;

    tmp_index[level] = 0;   update_node_HD(tmp_index, level+1, page_memory);
    tmp_index[level] = 1;   update_node_HD(tmp_index, level+1, page_memory);
    tmp_index[level] = 2;   update_node_HD(tmp_index, level+1, page_memory);
}

double KernalManager::bhv_PK427(std::vector<unsigned int> index, unsigned int level){
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    double Z0, Z1, Z2;

    tmp_index[level] = 0;   Z0 = get_bhattacharyya_value(tmp_index, level+1);
    tmp_index[level] = 1;   Z1 = get_bhattacharyya_value(tmp_index, level+1);
    tmp_index[level] = 2;   Z2 = get_bhattacharyya_value(tmp_index, level+1);
    
    Z0 = exp(-exp(Z0));
    Z1 = exp(-exp(Z1));
    Z2 = exp(-exp(Z2));

    switch(index[level]){
        case 0:
            // Z0+Z2-Z0Z2
            return log(-log(Z0+Z2-Z0*Z2));
            break;
        case 1:
            // Z1 + Z0Z2 - Z0Z1Z2
            return log(-log(Z1+Z0*Z2-Z0*Z1*Z2));
            break;
        case 2:
            // Z0Z1Z2
            return log(-log(Z0*Z1*Z2));
            break;
        default:
            ERROR("bhv_PK427 index with %d\n",index[level]);
            break;
    }   
    return 0;
}

void KernalManager::dnv_PK657(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1, site_2;
    double sub_node = 0;

    switch(index[level]){
        case 0:
            // need to get source
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);  do_node_value(tmp_index, level+1, page_memory);

            // u0 = x0x2
            page_memory[level][site_0].value = SIGN(page_memory[level+1][site_0].value) 
                                             * SIGN(page_memory[level+1][site_2].value)
                                             * MIN(ABS(page_memory[level+1][site_0].value),
                                                   ABS(page_memory[level+1][site_2].value));
            break;
        case 1:
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);

            // u1 = x1 (u0x2 or x0)
            sub_node = (page_memory[level][site_0].HD ? -1:1) * page_memory[level+1][site_2].value
                     + page_memory[level+1][site_0].value;

            page_memory[level][site_1].value = SIGN(page_memory[level+1][site_1].value) 
                                             * SIGN(sub_node)
                                             * MIN(ABS(page_memory[level+1][site_1].value),
                                                   ABS(sub_node));

            break;
        case 2:
            tmp_index[level] = 0;   site_0 = get_index(tmp_index);
            tmp_index[level] = 1;   site_1 = get_index(tmp_index);
            tmp_index[level] = 2;   site_2 = get_index(tmp_index);

            // u2 = u0u1x0 or u0x1 or  u1x2 
            page_memory[level][site_2].value = (page_memory[level][site_0].HD ? -1:1) * (page_memory[level][site_1].HD ? -1:1) * page_memory[level+1][site_0].value
                                             + (page_memory[level][site_0].HD ? -1:1) * page_memory[level+1][site_1].value
                                             + (page_memory[level][site_1].HD ? -1:1) * page_memory[level+1][site_2].value;

            break;
        default:
            ERROR("dnv_PK657 index with %d\n",index[level]);
            break;
    }
}

void KernalManager::unh_PK657(std::vector<unsigned int> index, unsigned int level, page &page_memory){
    // not ready update
    if(index[level] != 2)
        return; 

    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    unsigned int site_0, site_1, site_2;

    tmp_index[level] = 0;   site_0 = get_index(tmp_index);
    tmp_index[level] = 1;   site_1 = get_index(tmp_index);
    tmp_index[level] = 2;   site_2 = get_index(tmp_index);

    page_memory[level+1][site_0].HD = page_memory[level][site_0].HD 
                                    ^ page_memory[level][site_1].HD
                                    ^ page_memory[level][site_2].HD;
    page_memory[level+1][site_1].HD = page_memory[level][site_0].HD 
                                    ^ page_memory[level][site_2].HD;
    page_memory[level+1][site_2].HD = page_memory[level][site_1].HD 
                                    ^ page_memory[level][site_2].HD;

    tmp_index[level] = 0;   update_node_HD(tmp_index, level+1, page_memory);
    tmp_index[level] = 1;   update_node_HD(tmp_index, level+1, page_memory);
    tmp_index[level] = 2;   update_node_HD(tmp_index, level+1, page_memory);
}

double KernalManager::bhv_PK657(std::vector<unsigned int> index, unsigned int level){
    std::vector<unsigned int> tmp_index;
    tmp_index.assign(index.begin(),index.end());
    double Z0, Z1, Z2;

    tmp_index[level] = 0;   Z0 = get_bhattacharyya_value(tmp_index, level+1);
    tmp_index[level] = 1;   Z1 = get_bhattacharyya_value(tmp_index, level+1);
    tmp_index[level] = 2;   Z2 = get_bhattacharyya_value(tmp_index, level+1);
    
    Z0 = exp(-exp(Z0));
    Z1 = exp(-exp(Z1));
    Z2 = exp(-exp(Z2));

    switch(index[level]){
        case 0:
            // Z0+Z2-Z0Z2
            return log(-log(Z0+Z2-Z0*Z2));
            break;
        case 1:
            // Z1 + Z0Z2 - Z0Z1Z2
            return log(-log(Z1+Z0*Z2-Z0*Z1*Z2));
            break;
        case 2:
            // Z0Z1Z2
            return log(-log(Z0*Z1*Z2));
            break;
        default:
            ERROR("bhv_PK657 index with %d\n",index[level]);
            break;
    }   
    return 0;
}