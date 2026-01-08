#include "PolarMultiKernal_AED.h"
#include "BCH.h"
#include "libDebug.h"
#include "libParser.h"
#include "libMath.h"
#include <cmath>
#include <fstream>

PolarMultiKernal::PolarMultiKernal(std::map<std::string, std::string> config) : ErrorCorrectionCodeBase(config){
    // default setting
    ECC_info = "PolarMultiKernal";
    encoder_enable = false;
    decoder_enable = false;
    permutation_src = "";
    dynamic_frozen_process = "frozen";
    permutation_random_seed = -1;
    kernal_string = "";
    bha_value_setting = "";
    target_raw_BER = 0.1;
    list_size = 0;

    // AED
    useAED = false;
    aed_L = 0;
    automorphism_src = "";

    KM = NULL;
    stage = 0;

    received_order.assign(0,0);
    diverge_flag.assign(0,0);

    // deal config
    for(std::map<std::string, std::string>::iterator pair_idx = config.begin(); pair_idx != config.end(); pair_idx++){
        if      (pair_idx->first == "permutation_src")          permutation_src = pair_idx->second;
        else if (pair_idx->first == "dynamic_frozen_process")   dynamic_frozen_process = pair_idx->second;
        else if (pair_idx->first == "permutation_random_seed")  permutation_random_seed = stol(pair_idx->second);
        else if (pair_idx->first == "target_raw_BER")           target_raw_BER = stod(pair_idx->second);
        else if (pair_idx->first == "list_size")                list_size = stoi(pair_idx->second);
        else if (pair_idx->first == "kernal_string")            kernal_string = pair_idx->second;
        else if (pair_idx->first == "bha_value_setting")        bha_value_setting = pair_idx->second;
        else if (pair_idx->first == "useAED")                   useAED = (pair_idx->second == "1" || pair_idx->second == "true" || pair_idx->second == "True");
        else if (pair_idx->first == "use_AED")                  useAED = (pair_idx->second == "1" || pair_idx->second == "true" || pair_idx->second == "True");
        else if (pair_idx->first == "aed_L")                    aed_L = stoi(pair_idx->second);
        else if (pair_idx->first == "automorphism_src")         automorphism_src = pair_idx->second;
        else if (pair_idx->first == "message_length" && message_length == 0)            message_length = stoi(pair_idx->second);
        else {WARRING(" key : %s is not find in PolarMultiKernal.cpp.",pair_idx->first.c_str());}
    }
    // assert
    if(kernal_string == "")
        ERROR("kernal_string == \"\"");
    KM = new KernalManager(kernal_string);
    stage = KM->get_stage();
    codeword_length = KM->get_codeword_length();
    KM->build_Gmatrix(polar_matrix);

    if(DEBUG_MODE){
        printf("--polar_matrix--\n");
        printMatrix(polar_matrix);
        printf("--polar_matrix done-\n");
    }

    if(codeword_length < message_length)
        ERROR("codeword_length(%d) < message_length(%d)",codeword_length,message_length);
    if(list_size < 1)
        ERROR("list_size < 1");

    received_order.resize(codeword_length);
    diverge_flag.assign(codeword_length,0);
    relation_ship.resize(codeword_length);
    for(unsigned int idx=0; idx<codeword_length; idx++)
        relation_ship[idx].resize(0);

    // get permutation matrix
    if(permutation_src == "random"){
        permutationMatrix(permutation_matrix,codeword_length,&permutation_random_seed);
    }else if (permutation_src == "byBCH"){
        std::map<std::string, std::string> BCH_config;
        BCH_config["matrix_src"] = "byAlgorithm";
        BCH_config["codeword_length"] = std::to_string(codeword_length-1);
        BCH_config["message_length"] = std::to_string(message_length);
        BCH ECC_BCH = BCH(BCH_config);
        INFO("clone BCH_H");
        std::vector<std::vector<char> >BCH_Hmatrix = ECC_BCH.getHmatrix();
        Hmatrix.resize(BCH_Hmatrix.size()+1);
        for(unsigned int idx=0; idx<BCH_Hmatrix.size(); idx++){
            Hmatrix[idx].assign(BCH_Hmatrix[idx].size()+1,0);
            for(unsigned int jdx=0; jdx<BCH_Hmatrix[idx].size(); jdx++)
                Hmatrix[idx][jdx] = BCH_Hmatrix[idx][jdx];
            Hmatrix[idx][BCH_Hmatrix[idx].size()] = 1;
        }
        Hmatrix[BCH_Hmatrix.size()].assign(BCH_Hmatrix[0].size()+1,0);
        Hmatrix[BCH_Hmatrix.size()][BCH_Hmatrix[0].size()] = 1;
        INFO("clone BCH_G");
        std::vector<std::vector<char> >BCH_Gmatrix = ECC_BCH.getGmatrix();
        Gmatrix.resize(BCH_Gmatrix.size());
        for(unsigned int idx=0; idx<BCH_Gmatrix.size(); idx++){
            Gmatrix[idx].assign(BCH_Gmatrix[idx].size()+1,0);
            for(unsigned int jdx=0; jdx<BCH_Gmatrix[idx].size(); jdx++)
                Gmatrix[idx][jdx] = BCH_Gmatrix[idx][jdx];
            Gmatrix[idx][BCH_Gmatrix[idx].size()] = 1;
        }
        INFO("Gmatrix[%u][%u]",Gmatrix.size(), Gmatrix[0].size());
        std::vector<unsigned int> BCH_DtoA = ECC_BCH.getDtoA();
        permutation_matrix.resize(codeword_length);
        for(unsigned int idx=0; idx<codeword_length; idx++){
            permutation_matrix[idx].assign(codeword_length,0);
        }
        for(unsigned int idx=0; idx<codeword_length -1; idx++){
            permutation_matrix[codeword_length -1 - BCH_DtoA[idx]][idx] = 1;
        }
        permutation_matrix[codeword_length -1][codeword_length -1] = 1;
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

    // set received_order
    for(unsigned int idx=0; idx<codeword_length; idx++)
        for(unsigned int jdx=0; jdx<codeword_length; jdx++)
            if(permutation_matrix[idx][jdx])
                received_order[jdx] = idx;   // outside to inside

    if (bha_value_setting.size() == 0){

    }else if (bha_value_setting.size() == codeword_length){
        for(unsigned idx=0; idx<codeword_length; idx++)
            KM->bha_value_type[received_order[idx]] = bha_value_setting[idx];
    }else{
        ERROR("bha_value_setting error, size(%d) != codeword_length(%d)\n",bha_value_setting.size(), codeword_length);
    }

    KM->set_target_BER(target_raw_BER);

    if(DEBUG_MODE){
        printf("--permutation_matrix--\n");
        printMatrix(permutation_matrix);
        printf("--permutation_matrix done-\n");
    }

    if(matrix_src == "byAlgorithm"){
        // pure polar code, TODO
        // build Gmatrix, Hmatrix
        ERROR("TODO");
        // matrixMultiplication(polar_matrix, permutation_matrix, Gmatrix);
        // findNullSpace(Gmatrix, Hmatrix);
        // Gmatrix.resize(message_length);

        // // build diverge_flag, relation_ship
        // for(unsigned int idx=0, G_idx=0; idx<codeword_length; idx++){
        //     if(KM->bhattacharyya_order[idx]<message_length){
        //         // simple info
        //         diverge_flag[idx] = 1;
        //         Gmatrix[G_idx].assign(Hmatrix[idx].begin(), Hmatrix[idx].end());
        //         G_idx++;
        //     }else{
        //         // simple frozen
        //         relation_ship[idx].assign(1,idx);
        //     }
        //     INFO("bha %3d, %3d, %d", idx, KM->bhattacharyya_order[idx], diverge_flag[idx]);

        // }
        // encoder_enable = true;
        // findNullSpace(Gmatrix, Hmatrix);

    }else if(matrix_src == "byGmatrix" || matrix_src == "byHmatrix"){
        if(codeword_length != Hmatrix.size())
            ERROR("KM codeword_length(%d) != Hmatrix.size()(%d)",codeword_length,Hmatrix.size());

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

        // build diverge_flag, relation_ship
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
        for(unsigned int idx=0; idx<codeword_length; idx++){
            if(relation_ship[idx].size() == 0)   // info
                diverge_flag[idx] = 1;

            if(relation_ship[idx].size()>1 && dynamic_frozen_process == "info")  // dynamic
                diverge_flag[idx] = 1;
        }

        if(DEBUG_MODE){
            printf("--matrix & relation_ship--\n");
            for(unsigned idx=0; idx<polar_matrix.size(); idx++){
                for(unsigned jdx=0; jdx<polar_matrix[idx].size(); jdx++)
                    printf("%d ",polar_matrix[idx][jdx]);
                switch(relation_ship[idx].size()){
                    case 0:
                        printf("info");
                        printf("info          ");
                        break;
                    case 1:
                        printf("frozen        ");
                        break;
                    default:
                        printf("dynamic frozen");
                        break;
                }
                printf(" %f\n",exp(-exp(KM->bhattacharyya_value[idx])));
            }
            printf("--matrix & relation_ship done-\n");
        }
    }else
        ERROR("PolarMultiKernal not support matrix_src = %s",matrix_src.c_str());

    for(unsigned int idx=0; idx<codeword_length; idx++)
        printf("idx %2d : %lf\n",idx, exp(-exp(KM->bhattacharyya_value[idx])));

    double total_bhattacharyya_value = 0;
    for(unsigned int idx=0; idx<codeword_length; idx++){
        if(relation_ship[idx].size() == 0){
            total_bhattacharyya_value += exp(-exp(KM->bhattacharyya_value[idx]));
            printf("= %3d, %8f\n", idx, exp(-exp(KM->bhattacharyya_value[idx])));
        }
    }
    printf("get_basic_bhattacharyya_value %lf\n", exp(-exp(KM->get_basic_bhattacharyya_value())));
    printf("total_bhattacharyya_value %lf ",total_bhattacharyya_value );
    for (unsigned int idx=0; idx<codeword_length; idx++){
        printf("%d", relation_ship[idx].size() > 1 ? 2:relation_ship[idx].size());
    }printf("\n");
    // // if(total_bhattacharyya_value < exp(-exp(KM->get_basic_bhattacharyya_value()))*message_length){
    // FILE *out=fopen("result.txt","a+t");
    // fprintf(out, "%ld, %lf, ",temp_permutation_random_seed,total_bhattacharyya_value );
    // for (unsigned int idx=0; idx<codeword_length; idx++){
    //     fprintf(out, "%d", relation_ship[idx].size() > 1 ? 2:relation_ship[idx].size());
    // }fprintf(out, "\n");
    // fclose(out);
    // // }
    // exit(0);

    // mem
    list_enable.resize(list_size);
    path_metric.resize(list_size);
    extend_path_metric.resize(list_size*2);
    SCL_mem.resize(list_size);
    for(unsigned int SCL_idx=0; SCL_idx<list_size; SCL_idx++){
        SCL_mem[SCL_idx].resize(stage+1);
        for(unsigned int stage_idx=0; stage_idx<stage+1; stage_idx++)
            SCL_mem[SCL_idx][stage_idx].resize(codeword_length);
    }
    extend_message.resize(codeword_length);

    // ===== AED setup (v1-style): fixed constraints, vary only received->leaf mapping =====
    received_order_set.clear();
    if(useAED){
        if(automorphism_src == "")
            ERROR("useAED=true but automorphism_src is empty");
        loadAEDOrders(automorphism_src);
        buildReceivedOrderSet();
        if(aed_L == 0) aed_L = received_order_set.size();
        if(aed_L > received_order_set.size()) aed_L = received_order_set.size();
        INFO("AED enabled: trials=%u (available=%u)", aed_L, (unsigned)received_order_set.size());
    }else{
        received_order_set.push_back(received_order);
    }

    decoder_enable = true;
    // parseBinaryMatrix(Gmatrix_path,Gmatrix);
    // dumpHmatrix("PMK");
    // dumpGmatrix("PMK");
    if(DEBUG_MODE){
        INFO("done constructor");
        // exit(0);
    }
}


PolarMultiKernal::~PolarMultiKernal(){
    delete KM;
}

bool PolarMultiKernal::doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log){
    INFO("doDecode");

    // Ensure output sized
    if(decoded_word.size() != received.size())
        decoded_word.assign(received.size(), 0);

    // Single-shot (no AED)
    if(!useAED){
        double best_metric = 0;
        bool ok = decodeOnceWithOrder(received, received_order, decoded_word, best_metric, log);
        (void)ok;
        return true;
    }

    // AED: multiple trials, pick best valid (parity-check) else best metric
    const unsigned int L = (aed_L == 0) ? (unsigned int)received_order_set.size() : aed_L;
    std::vector<char> best_any(decoded_word.size(),0), best_valid(decoded_word.size(),0);
    double best_any_metric = 1e100;
    double best_valid_metric = 1e100;
    bool have_any = false;
    bool have_valid = false;

    for(unsigned int ell=0; ell<L; ell++){
        std::vector<char> cand(decoded_word.size(),0);
        double metric = 0;
        bool ok = decodeOnceWithOrder(received, received_order_set[ell], cand, metric, log);
        if(!ok) continue;
        have_any = true;
        if(metric < best_any_metric){
            best_any_metric = metric;
            best_any = cand;
        }
        if(parityCheckZero(cand)){
            if(!have_valid || metric < best_valid_metric){
                have_valid = true;
                best_valid_metric = metric;
                best_valid = cand;
            }
        }
    }

    if(have_valid){
        decoded_word = best_valid;
        INFO("AED selected VALID candidate metric=%e", best_valid_metric);
        return true;
    }
    if(have_any){
        decoded_word = best_any;
        INFO("AED no valid candidate; selected best metric=%e", best_any_metric);
        return false;
    }
    // fallback: should not happen
    INFO("AED found no candidate");
    return false;
}

void PolarMultiKernal::infoProcess(unsigned int decode_idx){
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

void PolarMultiKernal::frozenProcess(unsigned int decode_idx){
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

void PolarMultiKernal::clonePage(unsigned int src, unsigned int dst){
    INFO("clonePage %2d to %2d",src,dst);
    for(unsigned int level_idx=0; level_idx<SCL_mem[src].size(); level_idx++){
        for(unsigned int idx=0; idx<SCL_mem[src][level_idx].size(); idx++){
            SCL_mem[dst][level_idx][idx].value = SCL_mem[src][level_idx][idx].value;
            SCL_mem[dst][level_idx][idx].HD = SCL_mem[src][level_idx][idx].HD;
        }
    }
}

// ===== AED helpers =====

void PolarMultiKernal::loadAEDOrders(const std::string &path){
    aed_h_list.clear();
    if(path.size() == 0)
        return;

    std::ifstream fin(path.c_str());
    if(!fin.is_open())
        ERROR("cannot open automorphism_src = %s", path.c_str());

    // Supported formats:
    //  (1) header: <rows> <cols> followed by rows*cols integers, cols == N
    //  (2) no header: each line/row is N integers until EOF
    std::vector<unsigned int> buf;
    buf.reserve(codeword_length);

    // peek first two ints
    long a = -1, b = -1;
    if(!(fin >> a)) return;
    if(!(fin >> b)) ERROR("automorphism file truncated");

    bool has_header = false;
    const unsigned int rows = (unsigned int)a;
    const unsigned int cols = (unsigned int)b;
    if(cols == codeword_length){
        // Heuristic: treat as header if remaining count matches rows*cols OR
        // if rows seems to be a reasonable number of permutations.
        // We'll assume header.
        has_header = true;
    }

    if(has_header){
        for(unsigned int r=0; r<rows; r++){
            std::vector<unsigned int> h(codeword_length);
            for(unsigned int i=0; i<codeword_length; i++){
                unsigned int v;
                if(!(fin >> v)) ERROR("automorphism file truncated at row %u", r);
                if(v >= codeword_length) ERROR("automorphism value out of range: %u", v);
                h[i] = v;
            }
            aed_h_list.push_back(h);
        }
        return;
    }

    // no header: a,b were the first two entries of the first permutation row
    {
        std::vector<unsigned int> h(codeword_length);
        h[0] = (unsigned int)a;
        h[1] = (unsigned int)b;
        if(h[0] >= codeword_length || h[1] >= codeword_length)
            ERROR("automorphism value out of range in first row");
        for(unsigned int i=2; i<codeword_length; i++){
            unsigned int v;
            if(!(fin >> v)) ERROR("automorphism file truncated in first row");
            if(v >= codeword_length) ERROR("automorphism value out of range: %u", v);
            h[i] = v;
        }
        aed_h_list.push_back(h);
    }

    while(true){
        std::vector<unsigned int> h(codeword_length);
        for(unsigned int i=0; i<codeword_length; i++){
            unsigned int v;
            if(!(fin >> v)) return;
            if(v >= codeword_length) ERROR("automorphism value out of range: %u", v);
            h[i] = v;
        }
        aed_h_list.push_back(h);
    }
}

void PolarMultiKernal::buildReceivedOrderSet(){
    received_order_set.clear();
    // trial 0: baseline (no automorphism)
    received_order_set.push_back(received_order);

    for(unsigned int ell=0; ell<aed_h_list.size(); ell++){
        const std::vector<unsigned int> &h = aed_h_list[ell];
        if(h.size() != codeword_length) continue;
        std::vector<unsigned int> order(codeword_length);
        // h is NOT inverted: h[i] = h(i)
        for(unsigned int i=0; i<codeword_length; i++){
            order[i] = received_order[h[i]];
        }

        // de-duplicate
        bool dup = false;
        for(unsigned int k=0; k<received_order_set.size(); k++){
            if(received_order_set[k] == order){ dup = true; break; }
        }
        if(!dup) received_order_set.push_back(order);
    }
}

bool PolarMultiKernal::parityCheckZero(const std::vector<char> &codeword) const{
    if(Hmatrix.size() == 0)
        return true; // if no H provided, treat as pass

    // In this implementation we assume Hmatrix is N x (N-K)
    if(Hmatrix.size() != codeword_length)
        return false;
    if(Hmatrix[0].size() == 0)
        return true;

    for(unsigned int col=0; col<Hmatrix[0].size(); col++){
        char s = 0;
        for(unsigned int i=0; i<codeword_length; i++){
            s ^= (codeword[i] & Hmatrix[i][col]);
        }
        if(s) return false;
    }
    return true;
}

bool PolarMultiKernal::decodeOnceWithOrder(const std::vector<double> &received,
                                          const std::vector<unsigned int> &order,
                                          std::vector<char> &decoded_word,
                                          double &best_metric,
                                          std::string &log)
{
    (void)log;
    // reset list state
    list_enable.assign(list_size,0);
    path_metric.assign(list_size,0);
    extend_path_metric.assign(list_size*2,0);

    // clear SCL_mem (avoid stale state across AED trials)
    for(unsigned int SCL_idx=0; SCL_idx<list_size; SCL_idx++){
        for(unsigned int stage_idx=0; stage_idx<stage+1; stage_idx++){
            for(unsigned int i=0; i<codeword_length; i++){
                SCL_mem[SCL_idx][stage_idx][i].value = 0.0;
                SCL_mem[SCL_idx][stage_idx][i].HD = 0;
            }
        }
    }

    // load leaves using the provided outside->inside order map
    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++){
        const unsigned int leaf = order[codeword_idx];
        if (KM->bha_value_type[leaf] == "?")
            SCL_mem[0][stage][leaf].value = received[codeword_idx];
        else if (KM->bha_value_type[leaf] == "1")
            SCL_mem[0][stage][leaf].value = 99999999999999999.9999;
        else if (KM->bha_value_type[leaf] == "0")
            SCL_mem[0][stage][leaf].value = 0.0;
        else
            ERROR("bha_value_type[%d]=%s", codeword_idx, KM->bha_value_type[leaf]);
    }
    list_enable[0] = 1;

    // main SCL loop
    for(unsigned int decode_idx=0; decode_idx<codeword_length; decode_idx++){
        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx])
                KM->do_node_value(decode_idx, 0, SCL_mem[list_idx]);

        if(diverge_flag[decode_idx])
            infoProcess(decode_idx);
        else
            frozenProcess(decode_idx);

        for(unsigned int list_idx=0; list_idx<list_size; list_idx++)
            if(list_enable[list_idx])
                KM->update_node_HD(decode_idx, 0, SCL_mem[list_idx]);
    }

    // keep your original relation_ship-based pruning
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        for(unsigned int codeword_idx=0; codeword_idx<codeword_length && list_enable[list_idx]; codeword_idx++){
            if(diverge_flag[codeword_idx]){
                for(unsigned int relation_idx=0; relation_idx<relation_ship[codeword_idx].size(); relation_idx++)
                    list_enable[list_idx] ^= SCL_mem[list_idx][0][relation_ship[codeword_idx][relation_idx]].HD;
            }
        }
    }

    // pick min path_metric among surviving paths
    best_metric = 1e100;
    int target_list_idx = -1;
    for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
        if(list_enable[list_idx] && path_metric[list_idx] < best_metric){
            best_metric = path_metric[list_idx];
            target_list_idx = (int)list_idx;
        }
    }
    if(target_list_idx < 0)
        return false;

    if(decoded_word.size() != received.size())
        decoded_word.assign(received.size(),0);
    for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++)
        decoded_word[codeword_idx] = SCL_mem[target_list_idx][stage][order[codeword_idx]].HD;

    return true;
}
