#ifndef CHANNEL
#define CHANNEL
#include "struct.h"

#include <cstdint>
#include <vector>
#include "struct.h"
#include <random>

namespace PoAwN
{
    namespace channel
    {
        using structures::decoder_parameters;
        using structures::decoder_t;
        using structures::softdata_t;
        using structures::table_GF;
        using structures::vector;
        void EncodeChanBPSK_BinCCSK(std::mt19937 &gen, decoder_parameters &dec_param,
                                    const table_GF &table,
                                    const float SNR,
                                    const vector<vector<uint16_t>> &bin_table,
                                    vector<decoder_t> &chan_LLR_sorted,
                                    vector<uint16_t> &KSYMB,
                                    const vector<vector<softdata_t>> &bin_mod_dict);

        void EncodeChanGF_CCSK(std::mt19937 &gen, decoder_parameters &dec_param,
                               const table_GF &table,
                               const float SNR,
                               const vector<vector<uint16_t>> &CCSK_rotated_codes,
                               vector<decoder_t> &chan_LLR_sorted,
                               vector<uint16_t> &KSYMB);
    }
}

#endif