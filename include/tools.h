#ifndef TOOLS
#define TOOLS

#include "struct.h"
#include <cstdint>
#include <vector>

namespace PoAwN
{
    namespace tools
    {

        using std::vector;
        using namespace PoAwN::structures;

        void circshift(std::vector<uint16_t> &CCSK_seq, int shift);

        void create_ccsk_rotated_table(const vector<uint16_t> &CCSK_seq, const uint16_t q, vector<vector<uint16_t>> &ccsk_rotated_table);

        void RandomBinaryGenerator(const uint16_t K,
                                   const uint16_t q,
                                   const vector<vector<uint16_t>> &bin_table,
                                   const bool repeatable,
                                   const int SEED,
                                   vector<vector<uint16_t>> &KBIN,
                                   vector<uint16_t> &KSYMB);
        void RandomSymbGenerator(const uint16_t K,
                                 const uint16_t q,
                                 const bool repeatable,
                                 const int SEED,
                                 vector<uint16_t> &KSYMB);
        void AWGN_gen(double MEAN,
                      double STD,
                      bool repeatable,
                      double SEED,
                      vector<vector<softdata_t>> &noise_table);

        void awgn_channel_noise(const double sigma,
                                const bool repeatable,
                                const double SEED,
                                vector<vector<softdata_t>> &noisy_sig);
        void Encoder(const vector<vector<uint16_t>> &ADDGF, const vector<vector<uint16_t>> &MULGF,
                     const vector<vector<uint16_t>> &polar_coeff,
                     vector<vector<uint16_t>> &ucap, vector<uint16_t> &NSYMB);
        float My_drand48(int *initialise);
        // void inv_Encoder(const vector<vector<uint16_t>> &ADDGF, const vector<vector<uint16_t>> &DIVGF,
        //                                const vector<vector<uint16_t>> &polar_coeff,
        //                                const vector<uint16_t> NSYMB, vector<uint16_t> &u_symb);

    } // namespace tools
} // namespace PoAwN
#endif