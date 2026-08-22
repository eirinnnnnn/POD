#ifndef _MCLASS_DECODER_WRAPPER_H_
#define _MCLASS_DECODER_WRAPPER_H_

#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "AED_relation_check.h"
#include "libMath.h"

struct MClassTracePoint {
    unsigned int step = 0;
    unsigned int active = 0;
    double pm_min = 0.0;
    double pm_second = 0.0;
    double pm_gap = 0.0;
    double pm_mean = 0.0;
    double pm_max = 0.0;
    double llr_abs_min = 0.0;
    double llr_abs_mean = 0.0;
    double llr_abs_max = 0.0;
};

class MClassDecoder : public AdjustPolarDecoderRelation {
public:
    explicit MClassDecoder(std::map<std::string, std::string> config)
    : AdjustPolarDecoderRelation(config) {}

    unsigned int branchCount() const {
        return (unsigned int)received_order_set.size();
    }

    unsigned int messageLength() const {
        return message_length;
    }

    unsigned int codewordLength() const {
        return codeword_length;
    }

    unsigned int sclListSize() const {
        return list_size;
    }

    bool decodeBranch(const std::vector<double> &received,
                      unsigned int branch_idx,
                      std::vector<char> &candidate,
                      double &metric,
                      std::string &log) {
        if(branch_idx >= received_order_set.size())
            return false;
        return decode_with_order_scl(received,
                                     received_order_set[branch_idx],
                                     candidate,
                                     metric,
                                     log);
    }

    bool parityValid(std::vector<char> candidate) {
        std::vector<char> syndrome;
        matrixMultiplication(candidate, Hmatrix, syndrome);
        return oneCount(syndrome) == 0;
    }

    const std::vector<unsigned int>& branchOrder(unsigned int branch_idx) const {
        return received_order_set[branch_idx];
    }

    std::vector<unsigned int> relationSizes() const {
        std::vector<unsigned int> sizes(relation_ship.size(), 0);
        for(unsigned int i=0; i<relation_ship.size(); i++)
            sizes[i] = (unsigned int)relation_ship[i].size();
        return sizes;
    }

    // Raw dynamic-frozen relation row for decode_idx: the SCL_mem indices
    // (earlier decoded positions) whose hard decisions are XORed together to
    // force this position's value. Empty for a genuine information bit. By
    // construction the last entry equals decode_idx itself.
    const std::vector<unsigned int>& relationList(unsigned int decode_idx) const {
        return relation_ship[decode_idx];
    }

    bool decodeBranchTrace(const std::vector<double> &received,
                           unsigned int branch_idx,
                           const std::vector<unsigned int> &checkpoints,
                           std::vector<char> &candidate,
                           double &metric,
                           std::vector<MClassTracePoint> &trace) {
        trace.clear();
        if(branch_idx >= received_order_set.size())
            return false;
        const std::vector<unsigned int> &order = received_order_set[branch_idx];
        if(order.size() != codeword_length)
            return false;

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

        for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++){
            unsigned int target = order[codeword_idx];
            if(target >= codeword_length)
                return false;
            if (bha_value_setting[codeword_idx] == '?')
                SCL_mem[0][stage][target].value = received[codeword_idx];
            else if (bha_value_setting[codeword_idx] == '1')
                SCL_mem[0][stage][target].value = 99999999999999999.9999;
            else if (bha_value_setting[codeword_idx] == '0')
                SCL_mem[0][stage][target].value = 0.0;
            else
                return false;
        }

        std::vector<unsigned int> cps = checkpoints;
        std::sort(cps.begin(), cps.end());
        cps.erase(std::unique(cps.begin(), cps.end()), cps.end());
        unsigned int cp_idx = 0;

        auto collect = [&](unsigned int step) {
            MClassTracePoint p;
            p.step = step;
            double pm_sum = 0.0;
            double pm_min = std::numeric_limits<double>::max();
            double pm_second = std::numeric_limits<double>::max();
            double pm_max = -std::numeric_limits<double>::max();
            for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
                if(!list_enable[list_idx])
                    continue;
                p.active++;
                double pm = path_metric[list_idx];
                pm_sum += pm;
                if(pm < pm_min) {
                    pm_second = pm_min;
                    pm_min = pm;
                } else if(pm < pm_second) {
                    pm_second = pm;
                }
                if(pm > pm_max)
                    pm_max = pm;
            }
            if(p.active == 0) {
                p.pm_min = p.pm_second = p.pm_gap = p.pm_mean = p.pm_max = 0.0;
            } else {
                p.pm_min = pm_min;
                p.pm_second = (p.active >= 2 ? pm_second : pm_min);
                p.pm_gap = p.pm_second - p.pm_min;
                p.pm_mean = pm_sum / p.active;
                p.pm_max = pm_max;
            }

            unsigned int upto = std::min<unsigned int>(step, codeword_length);
            double abs_sum = 0.0;
            double abs_min = std::numeric_limits<double>::max();
            double abs_max = 0.0;
            unsigned int used = 0;
            for(unsigned int i=0; i<received.size(); i++) {
                if(order[i] >= upto)
                    continue;
                double v = std::abs(received[i]);
                abs_sum += v;
                if(v < abs_min)
                    abs_min = v;
                if(v > abs_max)
                    abs_max = v;
                used++;
            }
            if(used == 0) {
                p.llr_abs_min = p.llr_abs_mean = p.llr_abs_max = 0.0;
            } else {
                p.llr_abs_min = abs_min;
                p.llr_abs_mean = abs_sum / used;
                p.llr_abs_max = abs_max;
            }
            trace.push_back(p);
        };

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

            unsigned int step = decode_idx + 1;
            while(cp_idx < cps.size() && cps[cp_idx] == step) {
                collect(step);
                cp_idx++;
            }
        }
        while(cp_idx < cps.size()) {
            collect(cps[cp_idx]);
            cp_idx++;
        }

        for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
            for(unsigned int codeword_idx=0; codeword_idx<codeword_length && list_enable[list_idx]; codeword_idx++){
                if(diverge_flag[codeword_idx]){
                    for(unsigned int relation_idx=0; relation_idx<relation_ship[codeword_idx].size(); relation_idx++)
                        list_enable[list_idx] ^= SCL_mem[list_idx][0][relation_ship[codeword_idx][relation_idx]].HD;
                }
            }
        }

        double best = std::numeric_limits<double>::max();
        unsigned int target_list_idx = 0;
        bool any_enable = false;
        for(unsigned int list_idx=0; list_idx<list_size; list_idx++){
            if(list_enable[list_idx]){
                any_enable = true;
                if(path_metric[list_idx] < best){
                    best = path_metric[list_idx];
                    target_list_idx = list_idx;
                }
            }
        }
        if(!any_enable)
            return false;
        metric = path_metric[target_list_idx];
        candidate.resize(codeword_length);
        for(unsigned int codeword_idx=0; codeword_idx<received.size(); codeword_idx++)
            candidate[codeword_idx] = SCL_mem[target_list_idx][stage][order[codeword_idx]].HD;
        return true;
    }

    std::vector<unsigned int> informationPositions() const {
        std::vector<unsigned int> positions;
        for(unsigned int i=0; i<relation_ship.size(); i++) {
            if(relation_ship[i].empty())
                positions.push_back(i);
        }
        return positions;
    }
};

#endif
