#ifndef _POLAR_MULTI_KERNAL_PS_PED_H_
#define _POLAR_MULTI_KERNAL_PS_PED_H_
#include <map>
#include <string>
#include <vector>
#include "ErrorCorrectionCodeBase.h"

#include "KernalManager.h"

// Per-checkpoint snapshot of a branch's SCL list state (mirrors the eBCH
// side's MClassTracePoint). llr_abs_* are the raw channel-input magnitudes
// |received[i]| for codeword positions already "consumed" by this
// checkpoint (order[i] < checkpoint) -- decoder-agnostic, needs only
// `received` and the leaf-order map, not any internal SCL node state.
// Used to train/evaluate a branch-pruning classifier from a partial
// (mid-decode) signal instead of the final path metric.
struct PSPEDTracePoint {
    unsigned int active = 0;   // number of currently list_enable-alive paths
    double pm_min = 1e100;     // best (lowest) path_metric among alive paths
    double pm_gap = 0.0;       // second_best - best (0 if <2 alive paths)
    double pm_mean = 0.0;
    double pm_max = 0.0;
    double llr_abs_min = 0.0;
    double llr_abs_mean = 0.0;
    double llr_abs_max = 0.0;
};

class PolarMultiKernal : public ErrorCorrectionCodeBase{
public:
    PolarMultiKernal(std::map<std::string, std::string> config);
    virtual ~PolarMultiKernal();

    virtual bool doDecode(std::vector<double> &received, std::vector<char> &decoded_word, std::string &log);
    // virtual bool doEncode(std::vector<char> &message, std::vector<char> &codeword);
    
    // setTargetRBER will build and order bhattacharyya 
    // bool setTargetRBER(double RBER);

    // build relation_ship and made Gmatrix;
    // virtual bool buildRelationShip();


protected:
    // set by config
    std::string permutation_src;
    std::string dynamic_frozen_process; 
    long permutation_random_seed;
    std::string kernal_string;
    std::string bha_value_setting;
    double target_raw_BER;
    unsigned int list_size;

    KernalManager *KM;
    std::vector<std::vector<char> > polar_matrix;
    unsigned int stage;

    std::vector<std::vector<char> > permutation_matrix;
    std::vector<unsigned int> received_order;
    std::vector<char> diverge_flag;
    std::vector<std::vector<unsigned int> > relation_ship;
    // relation_ship[info bit idx].size = 0
    // relation_ship[(dynamic) forzen bit idx] last value is (dynamic) forzen bit idx

    // mem
    std::vector<char> list_enable;
    std::vector<double> path_metric;
    std::vector<double> extend_path_metric;
    std::vector<page> SCL_mem;
    std::vector<char> extend_message;

    // ===== AED (v1-style) =====
    // Toggle via config key: useAED=true/false (also accepts use_AED)
    // Automorphism file is an order-list (NOT inverted): h[i] = h(i)
    bool useAED;
    unsigned int aed_L;
    std::string automorphism_src;
    std::vector<std::vector<unsigned int>> aed_h_list;
    // composed order maps: order_ell[i] = received_order[h_ell[i]]
    std::vector<std::vector<unsigned int>> received_order_set;

    void loadAEDOrders(const std::string &path);
    void buildReceivedOrderSet();
    bool decodeOnceWithOrder(const std::vector<double> &received,
                             const std::vector<unsigned int> &order,
                             std::vector<char> &decoded_word,
                             double &best_metric,
                             std::string &log);
    bool parityCheckZero(const std::vector<char> &codeword) const;

public:
    unsigned int branchCount() const {
        if (!useAED) return 1;
        return aed_L == 0 ? (unsigned int)received_order_set.size() : aed_L;
    }

    // Decode every AED branch (up to L, or all available orders if L==0),
    // returning each branch's final candidate/metric/parity-validity PLUS a
    // running best-path-metric snapshot at each requested checkpoint
    // (checkpoints given as decode_idx counts, 1..codeword_length, ascending).
    // This is what makes genuine early-path-selection possible: ranking
    // branches by an EARLY checkpoint snapshot, then selecting among the
    // pruned survivors by their FINAL metric, is a lossy bet that can
    // diverge from the full-ensemble result -- ranking and selecting by the
    // same (final) metric, as decodeAllBranches() alone would do, cannot.
    bool decodeOnceWithOrderTrace(const std::vector<double> &received,
                                   const std::vector<unsigned int> &order,
                                   std::vector<char> &decoded_word,
                                   double &best_metric,
                                   std::string &log,
                                   const std::vector<unsigned int> &checkpoints,
                                   std::vector<PSPEDTracePoint> &checkpoint_trace);
    bool decodeAllBranchesTrace(const std::vector<double> &received,
                                 std::vector<std::vector<char> > &candidates,
                                 std::vector<double> &metrics,
                                 std::vector<char> &valid,
                                 const std::vector<unsigned int> &checkpoints,
                                 std::vector<std::vector<PSPEDTracePoint> > &checkpoint_trace,
                                 std::string &log);

    void infoProcess(unsigned int decode_idx);
    void frozenProcess(unsigned int decode_idx);
    void clonePage(unsigned int src, unsigned int dst);

    // decode_idx (0..codeword_length-1) -> true if this is a genuine info
    // bit (a free branching decision), false if frozen/dynamic-frozen
    // (deterministic given earlier decisions).
    const std::vector<char>& divergeFlags() const { return diverge_flag; }
    unsigned int relationSpan(unsigned int decode_idx) const { return (unsigned int)relation_ship[decode_idx].size(); }
    
};






#endif