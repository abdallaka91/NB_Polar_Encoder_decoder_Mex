#ifndef DECODER_FUNCTIONS
#define DECODER_FUNCTIONS
#include "struct.h"

#include <cstdint>
#include <vector>

namespace PoAwN
{

    namespace decoding
    {

        using structures::decoder_parameters;
        using structures::decoder_t;
        using structures::softdata_t;
        using structures::vector;

        void Channel_LLR(const vector<vector<softdata_t>> &chan_observ,
                         const vector<vector<uint16_t>> &bin_symb_seq,
                         uint16_t q,
                         softdata_t sigma,
                         vector<vector<softdata_t>> &chan_LLR);

        void LLR_sort(const vector<vector<softdata_t>> &chan_LLR,
                      const uint16_t nm,
                      vector<decoder_t> &chan_LLR_sorted,
                      const softdata_t ofst);

        // void ECN_EMS_bubble(const decoder_t &theta_1,
        //                     const decoder_t &phi_1,
        //                     const decoder_parameters &dec_param,
        //                     uint16_t coef,
        //                     const vector<vector<uint16_t>> &ADDGF,
        //                     const vector<vector<uint16_t>> &DIVGFC,
        //                     decoder_t &theta);
        void ECN_PA(const decoder_t &theta_1,
                    const decoder_t &phi_1,
                    const vector<vector<uint16_t>> &ADDGF,
                    const vector<vector<uint16_t>> &DIVGF,
                    const decoder_parameters &dec_param,
                    const uint16_t coef,
                    const uint16_t l,
                    const uint16_t s,
                    decoder_t &theta,
                    uint16_t &nm);

        void ECN_EMS(const decoder_t &theta_1,
                     const decoder_t &phi_1,
                     const vector<vector<uint16_t>> &ADDGF,
                     const vector<vector<uint16_t>> &DIVGF,
                     const decoder_parameters &dec_param,
                     const uint16_t coef,
                     decoder_t &theta);
        void ECN_EMS_L(const decoder_t &theta_1,
                       const decoder_t &phi_1,
                       const vector<vector<uint16_t>> &ADDGF,
                       const vector<vector<uint16_t>> &DIVGF,
                       const decoder_parameters &dec_param,
                       const uint16_t coef,
                       const uint16_t ucap_theta,
                       const uint16_t ucap_phi,
                       decoder_t &theta);

        void tuples_sorter(const decoder_t &tuples,
                           const uint16_t nm,
                           const uint16_t q,
                           const softdata_t offset,
                           decoder_t &sorted_tuples);

        void VN_update(const decoder_t &theta_1,
                       const decoder_t &phi_1,
                       const vector<vector<uint16_t>> &ADDGF,
                       const vector<vector<uint16_t>> &DIVGF,
                       const decoder_parameters &dec_param,
                       uint16_t coef,
                       uint16_t hard_decision,
                       decoder_t &phi);
        void VN_update_PA(const decoder_t &theta_1,
                          const decoder_t &phi_1,
                          const vector<vector<uint16_t>> &ADDGF,
                          const vector<vector<uint16_t>> &DIVGF,
                          const decoder_parameters &dec_param,
                          uint16_t coef,
                          const uint16_t l,
                          const uint16_t s,
                          const uint16_t df_nm,
                          uint16_t hard_decision,
                          decoder_t &phi);
        // void ECN_AEMS_bubble(const decoder_t &theta_1,
        //                      const decoder_t &phi_1,
        //                      const decoder_parameters &dec_param,
        //                      uint16_t coef,
        //                      const vector<vector<uint16_t>> &ADDGF,
        //                      const vector<vector<uint16_t>> &DIVGFC,
        //                      decoder_t &theta);
        void decode_SC_PA(const decoder_parameters &dec_param,
                          const vector<vector<uint16_t>> &ADDGF,
                          const vector<vector<uint16_t>> &MULGF,
                          const vector<vector<uint16_t>> &DIVGF,
                          vector<vector<decoder_t>> &L,
                          vector<uint16_t> &info_sec_rec);

        void decode_SC(const decoder_parameters &dec_param,
                                  const vector<vector<uint16_t>> &ADDGF,
                                  const vector<vector<uint16_t>> &MULGF,
                                  const vector<vector<uint16_t>> &DIVGF,
                                  vector<vector<decoder_t>> &L,
                                  vector<uint16_t> &info_sec_rec);

        void frozen_lay_pos(const decoder_parameters &dec_param,
                            vector<vector<uint16_t>> &ufrozen,
                            vector<vector<bool>> &clst_frozen);

    } // namespace decoding

} // namespace PoAwN

#endif