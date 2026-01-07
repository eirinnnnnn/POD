#include "AED_attemp_20260102.h"
#include "libDebug.h"
#include "libParser.h"
#include "libMath.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

AdjustPolarDecoder::AdjustPolarDecoder(std::map<std::string, std::string> config) : ErrorCorrectionCodeBase(config){
    // default setting
    ECC_info = "AdjustPolarDecoder";
    permutation_src = "";
    permutation_random_seed = -1;
    // dynamic_frozen_process = "frozen";
    target_raw_BER = 0.1;
    operationArray = "";
    list_size = 0;
    bha_value_setting = "";
    OnlyInit = false;

    // AED defaults
    use_AED = false;
    aed_L = 0;
    automorphism_src = "";

    stage = 0;

    received_order.assign(0,0);
    diverge_flag.assign(0,0);

    // deal config
    for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
        if      (pair_idx->first == "permutation_src")          permutation_src = pair_idx->second;
        else if (pair_idx->first == "permutation_random_seed")  permutation_random_seed = stol(pair_idx->second);
        // else if (pair_idx->first == "dynamic_frozen_process")   dynamic_frozen_process = pair_idx->second;
        else if (pair_idx->first == "target_raw_BER")           target_raw_BER = stod(pair_idx->second);
        else if (pair_idx->first == "operationArray")           operationArray = pair_idx->second;
        else if (pair_idx->first == "list_size")                list_size = stoi(pair_idx->second);
        else if (pair_idx->first == "bha_value_setting")        bha_value_setting = pair_idx->second;
        else if (pair_idx->first == "OnlyInit")                 OnlyInit = pair_idx->second == "true";
        else if (pair_idx->first == "use_AED")                  use_AED = (pair_idx->second == "true" || pair_idx->second == "True" || pair_idx->second == "1");
        else if (pair_idx->first == "aed_L")                    aed_L = stoi(pair_idx->second);
        else if (pair_idx->first == "automorphism_src")         automorphism_src = pair_idx->second;
        else {WARRING(" key : %s is not find in AdjustPolarDecoder.cpp.",pair_idx->first.c_str());}
    }
    // assert
    if(codeword_length & (codeword_length-1))
        ERROR("codeword_length(%d) is not power of 2", codeword_length);
    stage = uint(log2(codeword_length+0.5));
    if (operationArray.size() != stage*(codeword_length>>1))
        ERROR("operationArray error, size(%ld) != stage(%d)*codeword_length>>1(%d)\n", operationArray.size(), stage, (codeword_length>>1));
    if(list_size < 1)
        ERROR("list_size < 1");
    if (bha_value_setting.size() != codeword_length)
        ERROR("bha_value_setting error, size(%ld) != codeword_length(%d)\n",bha_value_setting.size(), codeword_length);

    // initial
    received_order.resize(codeword_length);
    diverge_flag.assign(codeword_length,0);
    relation_ship.resize(codeword_length);
    for(unsigned int idx=0; idx<codeword_length; idx++)
        relation_ship[idx].resize(0);
    list_enable.resize(list_size);
    path_metric.resize(list_size);
    extend_path_metric.resize(list_size*2);
    SCL_mem.resize(list_size);
    for(unsigned int SCL_idx=0; SCL_idx<list_size; SCL_idx++){
        SCL_mem[SCL_idx].resize(stage+1);
        for(unsigned int stage_idx=0; stage_idx<stage+1; stage_idx++)
            SCL_mem[SCL_idx][stage_idx].resize(codeword_length);
    }
    basic_bhattacharyya_value = log(-log(4*target_raw_BER*(1-target_raw_BER))/2);

    // get permutation matrix
    if(permutation_src == "random"){
        permutationMatrix(permutation_matrix,codeword_length,&permutation_random_seed);
    }else if(permutation_src != ""){
        parseBinaryMatrix(permutation_src,permutation_matrix);
    }else{
        printf("permutation matrix is idenity\n");
        permutation_matrix.resize(codeword_length);
        for(unsigned int idx=0; idx<codeword_length; idx++){
            permutation_matrix[idx].assign(codeword_length,0);
            permutation_matrix[idx][idx] = 1;
        }
    }
    if(DEBUG_MODE){
        printf("--permutation_matrix--\n");
        printMatrix(permutation_matrix);
        printf("--permutation_matrix done-\n");
    }

    // set received_order
    for(unsigned int idx=0; idx<codeword_length; idx++)
        for(unsigned int jdx=0; jdx<codeword_length; jdx++)
            if(permutation_matrix[idx][jdx])
                received_order[jdx] = idx;   // outside to inside

    // contruct polar_matrix
    polar_matrix.resize(codeword_length);
    for(unsigned int idx=0; idx<codeword_length; idx++){
        polar_matrix[idx].resize(codeword_length);
        polar_matrix[idx][idx] = 1;
    }

    for(unsigned int stage_idx=0; stage_idx<stage; stage_idx++){
        unsigned int idx=0, op_idx;
        unsigned int site_0, site_1;
        for(int jdx=0; jdx<int(codeword_length); jdx++){
            if(jdx & 0x1<<stage_idx){
                op_idx = stage_idx*(codeword_length>>1) + idx;
                idx++;
                site_0 = jdx ^ 0x1<<stage_idx;
                site_1 = jdx;
                if (operationArray[op_idx] == '0'){
                    for(unsigned int SCL_idx=0; SCL_idx<list_size; SCL_idx++){
                        SCL_mem[SCL_idx][stage_idx][site_0].linked = 0;
                        SCL_mem[SCL_idx][stage_idx][site_1].linked = 0;
                    }
                }else{
                    for(unsigned int SCL_idx=0; SCL_idx<list_size; SCL_idx++){
                        SCL_mem[SCL_idx][stage_idx][site_0].linked = 1;
                        SCL_mem[SCL_idx][stage_idx][site_1].linked = 1;
                    }
                    for(unsigned int tdx=site_1; tdx<codeword_length; tdx++){
                        if(polar_matrix[tdx][site_1])
                            polar_matrix[tdx][site_0] = 1;
                    }
                }
            }
        }
    }

    // bhattacharyya_value
    // in log(-log(Z)) domain, log(-log(4*target_BER*(1-target_BER))/2);
    bhattacharyya_value.resize(codeword_length);
    transmitted_length = 0;
    for(unsigned int idx=0; idx<codeword_length; idx++){
        if (bha_value_setting[idx] == '?'){
            transmitted_length++;
            bhattacharyya_value[received_order[idx]] = basic_bhattacharyya_value;
        }else if (bha_value_setting[idx] == '1')
            bhattacharyya_value[received_order[idx]] = 99999999999999999.9999;
        else if (bha_value_setting[idx] == '0')
            bhattacharyya_value[received_order[idx]] = log(-log(4*0.5*(1-0.5))/2);
        else
            ERROR("bha_value_setting[%d]=%c", idx, bha_value_setting[idx]);
    }
    for(int stage_idx=stage-1; stage_idx>=0; stage_idx--){
        unsigned int idx=0, op_idx;
        unsigned int site_0, site_1;
        for(int jdx=0; jdx<int(codeword_length); jdx++){
            if(jdx & 0x1<<stage_idx){
                op_idx = stage_idx*(codeword_length>>1) + idx;
                idx++;
                site_0 = jdx ^ 0x1<<stage_idx;
                site_1 = jdx;
                if (operationArray[op_idx] == '0'){
                    continue;
                }else{
                    double Z0 = exp(-exp(bhattacharyya_value[site_0]));
                    double Z1 = exp(-exp(bhattacharyya_value[site_1]));
                    bhattacharyya_value[site_0] = log(-log(Z0+Z1-Z0*Z1));
                    bhattacharyya_value[site_1] = log(-log(Z0*Z1));
                }
            }
        }
    }

    if(DEBUG_MODE){
        printf("--polar_matrix--\n");
        printMatrix(polar_matrix);
        printf("--polar_matrix linked & bhattacharyya_value--\n");
        for(int idx=0; idx<int(codeword_length); idx++){
            for(unsigned int stage_idx=0; stage_idx<stage; stage_idx++)
                printf("%d ", SCL_mem[0][stage_idx][idx].linked);
            printf("   -- %8f\n", exp(-exp(bhattacharyya_value[idx])));
        }
        printf("--polar_matrix done-\n");
    }


    if(matrix_src == "byGmatrix" || matrix_src == "byHmatrix"){
        if(DEBUG_MODE){
            printf("--Gmatrix--\n");
            printMatrix(Gmatrix);
            printf("--Gmatrix done-\n");
        }

        // cover Gmatrix
        std::vector<std::vector<char> > temp_Gmatrix, constraint, constraint_T;
        matrixMultiplication(polar_matrix, permutation_matrix, temp_Gmatrix);
        matrixMultiplication(temp_Gmatrix,Hmatrix,constraint);

        if(DEBUG_MODE){
            printf("--constraint--\n");
            printMatrix(constraint);
            printf("--constraint done-\n");
        }

        transposeMatrix(constraint, constraint_T);
        GaussianJordanElimination(constraint_T,1);

        if(DEBUG_MODE){
            printf("--GE constraint_T--\n");
            printMatrix(constraint_T);
            printf("--GE constraint_T done-\n");
        }

        // build relation_ship
        for(unsigned int row_idx=0; row_idx<constraint_T.size(); row_idx++){
            std::vector<unsigned int> tempV(0);
            for(unsigned int col_idx=0; col_idx<constraint_T[row_idx].size(); col_idx++){
                if(constraint_T[row_idx][col_idx])
                    tempV.push_back(col_idx);
            }
            if(tempV.size() == 0)
                continue;
            relation_ship[tempV[tempV.size()-1]].assign(tempV.begin(), tempV.end());
        }

        // build diverge_flag
        for(unsigned int idx=0; idx<codeword_length; idx++){
            if(relation_ship[idx].size() == 0)   // info
                diverge_flag[idx] = 1;
        }

        if(DEBUG_MODE){
            printf("--matrix & relation_ship--\n");
            for(unsigned idx=0; idx<polar_matrix.size(); idx++){
                for(unsigned jdx=0; jdx<polar_matrix[idx].size(); jdx++)
                    printf("%d ",polar_matrix[idx][jdx]);
                switch(relation_ship[idx].size()){
                    case 0:
                        printf("info          ");
                        break;
                    case 1:
                        printf("frozen        ");
                        break;
                    default:
                        printf("dynamic frozen");
                        break;
                }
                printf(" %f\n",exp(-exp(bhattacharyya_value[idx])));
            }
            printf("--matrix & relation_ship done-\n");
        }
    }else{
        ERROR("AdjustPolarDecoder not support matrix_src = %s",matrix_src.c_str());
    }

    // Prepare automorphism set for AED mode (always keep baseline order at index 0).
    if(!load_automorphism_set()){
        WARRING("load_automorphism_set failed, disable AED.");
        use_AED = false;
    }

    std::cout << "useAED" << use_AED << std::endl;
    std::cout << "received_order_set.empty()" << received_order_set.empty() << std::endl;
    for(unsigned int idx=0; idx<codeword_length; idx++)
        printf("idx %2d : %lf\n",idx, exp(-exp(bhattacharyya_value[idx])));

    double total_bhattacharyya_value = 0;
    for(unsigned int idx=0; idx<codeword_length; idx++){
        if(relation_ship[idx].size() == 0){
            total_bhattacharyya_value += exp(-exp(bhattacharyya_value[idx]));
            printf("= %3d, %8f\n", idx, exp(-exp(bhattacharyya_value[idx])));
        }
    }

    printf("get_basic_bhattacharyya_value %lf\n", exp(-exp(basic_bhattacharyya_value)));
    printf("total_bhattacharyya_value %lf ",total_bhattacharyya_value );
    for (unsigned int idx=0; idx<codeword_length; idx++){
        printf("%ld", relation_ship[idx].size() > 1 ? 2:relation_ship[idx].size());
    }printf("\n");

    if(OnlyInit)
        exit(0);
    INFO("done constructor");
}


AdjustPolarDecoder::~AdjustPolarDecoder(){

}

double AdjustPolarDecoder::getCodeRate(){
    return 1.0*message_length/transmitted_length;
}

bool AdjustPolarDecoder::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
    INFO("doDecode");

    if(use_AED){
        return doDecode_AED(received, decoded_word, log);
    }

    list_enable.assign(list_size,0);
    path_metric.assign(list_size,0);

    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++){
        if (bha_value_setting[codeword_idx] == '?')
            SCL_mem[0][stage][received_order[codeword_idx]].value = received[codeword_idx]; // need to change to LLR domain
        else if (bha_value_setting[codeword_idx] == '1')
            SCL_mem[0][stage][received_order[codeword_idx]].value = 99999999999999999.9999;
        else if (bha_value_setting[codeword_idx] == '0')
            SCL_mem[0][stage][received_order[codeword_idx]].value = 0.0;
        else
            ERROR("bha_value_setting[%d]=%c", codeword_idx, bha_value_setting[codeword_idx]);
    }
    list_enable[0] = 1;

    for(unsigned int decode_idx=0; decode_idx<codeword_length; decode_idx++){
        INFO("decode_idx %2d",decode_idx);
        // compute node
        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx]){
                do_node_value(decode_idx, 0, SCL_mem[list_idx]);
                INFO("SCL_mem[%2d][0][%2d] = %lf",list_idx,decode_idx,SCL_mem[list_idx][0][decode_idx].value);
            }
        INFO("");

        // diverge or not
        if(diverge_flag[decode_idx])
            infoProcess(decode_idx);
        else
            frozenProcess(decode_idx);

        // update node
        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx]){
                update_node_HD(decode_idx, 0, SCL_mem[list_idx]);
            }

        for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
            INFO("path_metric[%2d] = %e", list_idx, path_metric[list_idx]);
        }

        INFO("");
    }

    // check releation ship
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        for(unsigned int codeword_idx=0; codeword_idx<codeword_length && list_enable[list_idx]; codeword_idx++){
            if(diverge_flag[codeword_idx]){
                for(unsigned int relation_idx=0; relation_idx<relation_ship[codeword_idx].size(); relation_idx++)
                    list_enable[list_idx] ^= SCL_mem[list_idx][0][relation_ship[codeword_idx][relation_idx]].HD;
            }
        }
    }

    // pick min path_metric
    double tmp_double = 99999999;
    unsigned int target_list_idx = 0;
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        INFO("path_metric[%2d] = %e", list_idx, path_metric[list_idx]);
        if(list_enable[list_idx] && path_metric[list_idx]<tmp_double){
            target_list_idx = list_idx;
            tmp_double = path_metric[list_idx];
        }
    }
    // printf("path_metric %lf \n",path_metric[target_list_idx] );
    // deal output
    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++)
        decoded_word[codeword_idx] = SCL_mem[target_list_idx][stage][received_order[codeword_idx]].HD;

    return true;
}

bool AdjustPolarDecoder::doDecode_AED(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
    // std::cout << "decode AED\n";
    if(received_order_set.empty()){
        WARRING("AED enabled but no permutation set is loaded. Fallback to standard SCL.");
        use_AED = false;
        return doDecode(received, decoded_word, log);
    }

    unsigned int perm_limit = received_order_set.size();

    INFO("doDecode_AED: try %u permutations (base + automorphisms)", perm_limit);

    bool found_valid = false;
    double best_metric = std::numeric_limits<double>::max();
    std::vector<char> best_word(codeword_length, 0);

    for(unsigned int perm_idx=0; perm_idx<perm_limit; perm_idx++){
        std::vector<char> candidate(codeword_length, 0);
        double metric = 0.0;
        if(!decode_with_order_scl(received, received_order_set[perm_idx], candidate, metric, log))
            continue;

        std::vector<char> syndrome;
        matrixMultiplication(candidate, Hmatrix, syndrome);
        bool valid = (oneCount(syndrome) == 0);

        if(valid && (!found_valid || metric < best_metric)){
            best_metric = metric;
            best_word.swap(candidate);
            found_valid = true;
            // std::cout << "a working permutation\n";
        }else if(!found_valid && metric < best_metric){
            best_metric = metric;
            best_word.swap(candidate);
            // std::cout << "a working permutation\n";
        }
    }

    if(found_valid){
        decoded_word.swap(best_word);
        return true;
    }

    // No syndrome-clean candidate; return the most reliable one we saw.
    decoded_word.swap(best_word);
    return false;
}

bool AdjustPolarDecoder::decode_with_order_scl(const std::vector<double> &received,
                                               const std::vector<unsigned int> &order,
                                               std::vector<char> &decoded_word,
                                               double &metric,
                                               std::string &log){
    metric = 0.0;
    if(order.size() != codeword_length){
        WARRING("decode_with_order_scl: order size (%lu) != codeword_length (%u)", order.size(), codeword_length);
        decoded_word.assign(codeword_length, 0);
        return false;
    }

    // Reset working memory for all list entries.
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        for(unsigned int level_idx=0; level_idx<SCL_mem[list_idx].size(); level_idx++){
            for(unsigned int idx=0; idx<SCL_mem[list_idx][level_idx].size(); idx++){
                SCL_mem[list_idx][level_idx][idx].HD = 0;
                SCL_mem[list_idx][level_idx][idx].value = 0.0;
            }
        }
    }
    std::fill(list_enable.begin(), list_enable.end(), 0);
    std::fill(path_metric.begin(), path_metric.end(), 0.0);
    list_enable[0] = 1;

    // Apply permutation order to received LLRs.
    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++){
        unsigned int target = order[codeword_idx];
        if(target >= codeword_length){
            WARRING("decode_with_order_scl: target index %u out of range", target);
            return false;
        }
        if (bha_value_setting[codeword_idx] == '?')
            SCL_mem[0][stage][target].value = received[codeword_idx];
        else if (bha_value_setting[codeword_idx] == '1')
            SCL_mem[0][stage][target].value = 99999999999999999.9999;
        else if (bha_value_setting[codeword_idx] == '0')
            SCL_mem[0][stage][target].value = 0.0;
        else
            ERROR("bha_value_setting[%d]=%c", codeword_idx, bha_value_setting[codeword_idx]);
    }

    // SCL over the permuted order.
    for(unsigned int decode_idx=0; decode_idx<codeword_length; decode_idx++){
        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx])
                do_node_value(decode_idx, 0, SCL_mem[list_idx]);

        if(diverge_flag[decode_idx])
            infoProcess(decode_idx);
        else
            frozenProcess(decode_idx);

        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx])
                update_node_HD(decode_idx, 0, SCL_mem[list_idx]);
    }

    // Check dynamic frozen relation.
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        for(unsigned int codeword_idx=0; codeword_idx<codeword_length && list_enable[list_idx]; codeword_idx++){
            if(diverge_flag[codeword_idx]){
                for(unsigned int relation_idx=0; relation_idx<relation_ship[codeword_idx].size(); relation_idx++)
                    list_enable[list_idx] ^= SCL_mem[list_idx][0][relation_ship[codeword_idx][relation_idx]].HD;
            }
        }
    }

    // Pick best path.
    double tmp_double = std::numeric_limits<double>::max();
    unsigned int target_list_idx = 0;
    bool any_enable = false;
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        if(list_enable[list_idx]){
            any_enable = true;
            if(path_metric[list_idx] < tmp_double){
                tmp_double = path_metric[list_idx];
                target_list_idx = list_idx;
            }
        }
    }
    if(!any_enable){
        WARRING("decode_with_order_scl: all paths disabled after relation check.");
        return false;
    }
    metric = path_metric[target_list_idx];

    decoded_word.resize(codeword_length);
    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++)
        decoded_word[codeword_idx] = SCL_mem[target_list_idx][stage][order[codeword_idx]].HD;
    return true;
}

bool AdjustPolarDecoder::parse_order_from_block(const std::vector<std::vector<unsigned int> > &raw,
                                                unsigned int offset,
                                                std::vector<unsigned int> &order){
    order.assign(codeword_length, codeword_length);
    std::vector<char> col_seen(codeword_length, 0);

    for(unsigned int row=0; row<codeword_length; row++){
        unsigned int row_one_count = 0;
        unsigned int one_col = 0;
        for(unsigned int col=0; col<codeword_length; col++){
            unsigned int val = raw[offset + row][col];
            if(val == 1){
                row_one_count++;
                one_col = col;
            }else if(val != 0){
                return false;   // not a permutation matrix row
            }
        }
        if(row_one_count != 1 || col_seen[one_col])
            return false;
        col_seen[one_col] = 1;
        order[one_col] = row;
    }

    for(unsigned int idx=0; idx<order.size(); idx++)
        if(order[idx] == codeword_length)
            return false;
    return true;
}
// bool AdjustPolarDecoder::load_automorphism_set(){
//     received_order_set.clear();

//     if(!use_AED && automorphism_src == "")
//         return true;
//     if(automorphism_src == ""){
//         INFO("AED enabled without automorphism_src, using only baseline permutation.");
//         received_order_set.push_back(received_order);
//         return true;
//     }

//     FILE *fp = fopen(automorphism_src.c_str(), "r");
//     if(!fp){
//         WARRING("automorphism_src %s cannot be opened.", automorphism_src.c_str());
//         return false;
//     }
//     unsigned int rows = 0, cols = 0;
//     if(fscanf(fp, "%u %u", &rows, &cols) != 2){
//         WARRING("automorphism_src %s header format error.", automorphism_src.c_str());
//         fclose(fp);
//         return false;
//     }
//     std::vector<std::vector<unsigned int> > raw(rows, std::vector<unsigned int>(cols, 0));
//     for(unsigned int r=0; r<rows; r++){
//         for(unsigned int c=0; c<cols; c++){
//             if(fscanf(fp, "%u", &raw[r][c]) != 1){
//                 WARRING("automorphism_src %s body format error at row %u col %u.", automorphism_src.c_str(), r, c);
//                 fclose(fp);
//                 return false;
//             }
//         }
//     }
//     fclose(fp);

//     if(cols != codeword_length){
//         WARRING("automorphism_src columns (%u) != codeword_length (%u).", cols, codeword_length);
//         return false;
//     }

//     std::vector<std::vector<unsigned int> > orders;

//     // Try to parse as direct order list (each row is a permutation array).
//     bool as_order = true;
//     for(unsigned int r=0; r<rows; r++){
//         std::vector<char> used(codeword_length, 0);
//         for(unsigned int c=0; c<cols; c++){
//             unsigned int val = raw[r][c];
//             if(val >= codeword_length || used[val]){
//                 as_order = false;
//                 break;
//             }
//             used[val] = 1;
//         }
//         if(!as_order)
//             break;
//     }
//     if(as_order){
//         orders = raw;
//     }else if(rows % codeword_length == 0){
//         // Fallback: interpret as stacked permutation matrices.
//         unsigned int perm_cnt = rows / codeword_length;
//         for(unsigned int p=0; p<perm_cnt; p++){
//             std::vector<unsigned int> order;
//             if(parse_order_from_block(raw, p*codeword_length, order))
//                 orders.push_back(order);
//         }
//     }

//     if(orders.empty()){
//         WARRING("automorphism_src %s cannot be parsed into permutations.", automorphism_src.c_str());
//         return false;
//     }

//     // Compute inverse of received_order (base permutation)
//     // received_order[j] = i means position j receives from position i
//     // inverse_received_order[i] = j means position i sends to position j
//     std::vector<unsigned int> inverse_received_order(codeword_length, 0);
//     for(unsigned int j=0; j<codeword_length; j++){
//         inverse_received_order[received_order[j]] = j;
//     }

//     unsigned int limit = aed_L ? std::min<unsigned int>(aed_L, orders.size()) : orders.size();
//     for(unsigned int idx=0; idx<limit; idx++){
//         std::vector<unsigned int> composed(codeword_length, 0);
        
//         // New composition: h^{-1} ∘ P ∘ P^{-1} where P^{-1} is received_order
//         // Step 1: Apply inverse of base permutation (unpermute)
//         // Step 2: Apply identity (base permutation cancels out, so we just apply automorphism)
//         // Step 3: Apply automorphism
//         // Result: For each output position j, trace back through automorphism only
        
//         // for(unsigned int j=0; j<codeword_length; j++){
//         //     // Apply automorphism to the unpermuted position
//         //     unsigned int auto_input = inverse_received_order[j];  // unpermute position j
//         //     unsigned int auto_output = orders[idx][auto_input];   // apply automorphism
//         //     composed[j] = received_order[auto_output];             // re-permute
//         // }
//         for(unsigned int j=0; j<codeword_length; j++){
//             composed[j] = orders[idx][received_order[j]];  // h^{-1}[P^{-1}[j]]
// }

//         bool dup = false;
//         for(const auto &existing : received_order_set){
//             if(existing == composed){
//                 dup = true;
//                 break;
//             }
//         }
//         if(!dup)
//             received_order_set.push_back(composed);
//     }

//     if(received_order_set.empty()){
//         // Fallback to baseline permutation if nothing parsed.
//         received_order_set.push_back(received_order);
//     }

//     INFO("load_automorphism_set: using %lu permutations (including baseline).", (unsigned long)received_order_set.size());
//     return true;
// }
bool AdjustPolarDecoder::load_automorphism_set(){
    received_order_set.clear();

    if(!use_AED && automorphism_src == "")
        return true;
    if(automorphism_src == ""){
        INFO("AED enabled without automorphism_src, using only baseline permutation.");
        received_order_set.push_back(received_order);
        return true;
    }

    FILE *fp = fopen(automorphism_src.c_str(), "r");
    if(!fp){
        WARRING("automorphism_src %s cannot be opened.", automorphism_src.c_str());
        return false;
    }
    unsigned int rows = 0, cols = 0;
    if(fscanf(fp, "%u %u", &rows, &cols) != 2){
        WARRING("automorphism_src %s header format error.", automorphism_src.c_str());
        fclose(fp);
        return false;
    }
    std::vector<std::vector<unsigned int> > raw(rows, std::vector<unsigned int>(cols, 0));
    for(unsigned int r=0; r<rows; r++){
        for(unsigned int c=0; c<cols; c++){
            if(fscanf(fp, "%u", &raw[r][c]) != 1){
                WARRING("automorphism_src %s body format error at row %u col %u.", automorphism_src.c_str(), r, c);
                fclose(fp);
                return false;
            }
        }
    }
    fclose(fp);

    if(cols != codeword_length){
        WARRING("automorphism_src columns (%u) != codeword_length (%u).", cols, codeword_length);
        return false;
    }

    std::vector<std::vector<unsigned int> > orders;

    // Try to parse as direct order list (each row is a permutation array).
    bool as_order = true;
    for(unsigned int r=0; r<rows; r++){
        std::vector<char> used(codeword_length, 0);
        for(unsigned int c=0; c<cols; c++){
            unsigned int val = raw[r][c];
            if(val >= codeword_length || used[val]){
                as_order = false;
                break;
            }
            used[val] = 1;
        }
        if(!as_order)
            break;
    }
    if(as_order){
        orders = raw;
    }else if(rows % codeword_length == 0){
        // Fallback: interpret as stacked permutation matrices.
        unsigned int perm_cnt = rows / codeword_length;
        for(unsigned int p=0; p<perm_cnt; p++){
            std::vector<unsigned int> order;
            if(parse_order_from_block(raw, p*codeword_length, order))
                orders.push_back(order);
        }
    }

    if(orders.empty()){
        WARRING("automorphism_src %s cannot be parsed into permutations.", automorphism_src.c_str());
        return false;
    }

    unsigned int limit = aed_L ? std::min<unsigned int>(aed_L, orders.size()) : orders.size();
    for(unsigned int idx=0; idx<limit; idx++){
        std::vector<unsigned int> composed(codeword_length, 0);
        for(unsigned int j=0; j<codeword_length; j++)
            composed[j] = received_order[orders[idx][j]];   // P^{-1} ∘ h^{-1}

        bool dup = false;
        for(const auto &existing : received_order_set){
            if(existing == composed){
                dup = true;
                break;
            }
        }
        if(!dup)
            received_order_set.push_back(composed);
    }

    if(received_order_set.empty()){
        // Fallback to baseline permutation if nothing parsed.
        received_order_set.push_back(received_order);
    }

    INFO("load_automorphism_set: using %lu permutations (including baseline).", (unsigned long)received_order_set.size());
    return true;
}

void AdjustPolarDecoder::do_node_value(unsigned int index, unsigned int level, page &page_memory){
    INFO("do_node_value level %d",level);
    if(level == stage)
        return;
    unsigned int site_0, site_1;
    if(index & 0x1<<level){
        // site 1
        site_0 = index ^ 0x1<<level;
        site_1 = index;
        if(page_memory[level][site_1].linked){
            page_memory[level][site_1].value =  page_memory[level+1][site_1].value
                                             + (page_memory[level  ][site_0].HD ? -1:1)
                                             *  page_memory[level+1][site_0].value;
        }else{
            page_memory[level][site_1].value = page_memory[level+1][site_1].value;
        }
    }else{
        // site 0
        site_0 = index;
        site_1 = index ^ 0x1<<level;
        do_node_value(site_0, level+1, page_memory);
        do_node_value(site_1, level+1, page_memory);
        if(page_memory[level][site_0].linked){
            page_memory[level][site_0].value = SIGN(page_memory[level+1][site_0].value)
                                             * SIGN(page_memory[level+1][site_1].value)
                                             * MIN(ABS(page_memory[level+1][site_0].value),
                                                   ABS(page_memory[level+1][site_1].value));
        }else{
            page_memory[level][site_0].value = page_memory[level+1][site_0].value;
        }
    }
}

void AdjustPolarDecoder::update_node_HD(unsigned int index, unsigned int level, page &page_memory){
    INFO("update_node_HD level %d",level);
    if(level == stage)
        return;
    unsigned int site_0, site_1;
    if(index & 0x1<<level){
        // site 1
        site_0 = index ^ 0x1<<level;
        site_1 = index;
        if(page_memory[level][site_1].linked){
            page_memory[level+1][site_0].HD = page_memory[level][site_0].HD
                                            ^ page_memory[level][site_1].HD;
            page_memory[level+1][site_1].HD = page_memory[level][site_1].HD;
        }else{
            page_memory[level+1][site_0].HD = page_memory[level][site_0].HD;
            page_memory[level+1][site_1].HD = page_memory[level][site_1].HD;
        }
        update_node_HD(site_0, level+1, page_memory);
        update_node_HD(site_1, level+1, page_memory);
    }else{
        // site 0
        // not ready update
        return;
    }
}

void AdjustPolarDecoder::infoProcess(unsigned int decode_idx){
    extend_path_metric.assign(list_size*2,-1);
    unsigned int extend_list_count=0;

    // extend path
    for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
        if(list_enable[list_idx]){
            extend_list_count += 2;
            if(SCL_mem[list_idx][0][decode_idx].value < 0){
                extend_path_metric[list_idx*2  ] = path_metric[list_idx] + ABS(SCL_mem[list_idx][0][decode_idx].value);
                extend_path_metric[list_idx*2+1] = path_metric[list_idx];
            }else{
                extend_path_metric[list_idx*2  ] = path_metric[list_idx];
                extend_path_metric[list_idx*2+1] = path_metric[list_idx] + ABS(SCL_mem[list_idx][0][decode_idx].value);
            }
        }
    }

    // maintain path_metric.size() = list_size, kill large path_metric
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
                SCL_mem[list_idx][0][decode_idx].HD = 1;
            }else{
                path_metric[list_idx] = extend_path_metric[list_idx*2];
                SCL_mem[list_idx][0][decode_idx].HD = 0;
            }
        }
    }

}

void AdjustPolarDecoder::frozenProcess(unsigned int decode_idx){
    for(unsigned int list_idx=0; list_idx<list_enable.size(); list_idx++){
        if(list_enable[list_idx]){
            SCL_mem[list_idx][0][decode_idx].HD = 0;
            for(unsigned int relation_idx=0; relation_idx<relation_ship[decode_idx].size(); relation_idx++)
                if(decode_idx != relation_ship[decode_idx][relation_idx])
                    SCL_mem[list_idx][0][decode_idx].HD ^= SCL_mem[list_idx][0][relation_ship[decode_idx][relation_idx]].HD;
            if(SIGN(SCL_mem[list_idx][0][decode_idx].value) * (0.5- SCL_mem[list_idx][0][decode_idx].HD) < 0)
                path_metric[list_idx] += ABS(SCL_mem[list_idx][0][decode_idx].value);
        }
    }
}

void AdjustPolarDecoder::clonePage(unsigned int src, unsigned int dst){
    INFO("clonePage %2d to %2d",src,dst);
    for(unsigned int level_idx=0; level_idx<SCL_mem[src].size(); level_idx++){
        for(unsigned int idx=0; idx<SCL_mem[src][level_idx].size(); idx++){
            SCL_mem[dst][level_idx][idx].value = SCL_mem[src][level_idx][idx].value;
            SCL_mem[dst][level_idx][idx].HD = SCL_mem[src][level_idx][idx].HD;
        }
    }
}
