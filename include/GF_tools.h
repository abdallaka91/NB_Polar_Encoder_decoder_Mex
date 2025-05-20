#ifndef GF_TOOLS_H
#define GF_TOOLS_H

#include <vector>
#include <cstdint>

namespace PoAwN
{
    namespace GFtools
    {

        void GF_bin_seq_gen(uint16_t order, uint16_t primitivePoly,
                            std::vector<std::vector<uint16_t>> &GF_binSeq,
                            std::vector<std::vector<uint16_t>> &Dec_binSeq);
        void GF_bin2GF(const std::vector<std::vector<uint16_t>> &GF_bin,
                       std::vector<uint16_t> &dec_GF,
                       std::vector<uint16_t> &GF_dec);
        void GF_add_mat_gen(const std::vector<std::vector<uint16_t>> &BINGF,
                            const std::vector<uint16_t> &DECGF,
                            const std::vector<uint16_t> &GFDEC,
                            std::vector<std::vector<uint16_t>> &GF_add_mat,
                            std::vector<std::vector<uint16_t>> &DEC_add_mat);
        void GF_mul_mat_gen(const std::vector<uint16_t> &DECGF,
                            std::vector<std::vector<uint16_t>> &GF_mul_mat,
                            std::vector<std::vector<uint16_t>> &DEC_mul_mat);
        void GF_div_mat_gen(const std::vector<uint16_t> &DECGF,
                            std::vector<std::vector<uint16_t>> &GF_div_mat,
                            std::vector<std::vector<uint16_t>> &DEC_div_mat);
    } // namespace GFtools
}
#endif