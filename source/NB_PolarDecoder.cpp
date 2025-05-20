#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "Decoder_functions.h"
#include "GF_tools.h"
#include "init.h"
#include "struct.h"
#include "tools.h"
#include "HelperFunc.h"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>
#include <iomanip>
#include <algorithm>
#include "channel.h"
#include <filesystem>
#include <omp.h>
#include <atomic>
#include <random>
#include "mex.h"

// #include <omp.h>

using namespace PoAwN::structures;
using namespace PoAwN::tools;
using namespace PoAwN::init;
using namespace PoAwN::decoding;
using namespace PoAwN::channel;
using std::array;
using std::cout;
using std::endl;
using std::stod;
using std::stoi;
using std::string;
using std::vector;

void NBPolarDecoder(const std::vector<double> &LLR_in,
                    const uint16_t K,
                    const std::vector<uint16_t> &reliab_sequence,
                    const uint16_t nm,
                    const double ofst,
                    std::vector<uint16_t> &KSYMB)
{

    uint16_t N = reliab_sequence.size();
    KSYMB.resize(K, 0);
    string sig_mod = "CCSK_NB";
    uint16_t frozen_val = 0;
    uint16_t n = log2(N);
    uint16_t q = 64;
    uint16_t p = log2(q);
    softdata_t offset = (softdata_t)ofst;
    uint16_t nH = nm, nL = nm, Zc = 2, nb = 4, nopM = 4;

    base_code_t code_param(N, K, n, q, p, frozen_val);
    code_param.sig_mod = sig_mod;
    code_param.reliab_sequence = reliab_sequence;

    table_GF table;
    // LoadCode(code_param, EbN0);
    LoadTables(code_param, table, GF_polynom_primitive.data());
    code_param.polar_coeff.resize(n, std::vector<uint16_t>(N / 2));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < N / 2; j++)
        {
            code_param.polar_coeff[i][j] = 1;
        }
    }

    decoder_parameters dec_param(code_param, offset, nm, nL, nH, nb, Zc, nopM);

    dec_param.Roots_V.resize(n + 1);
    dec_param.Roots_indices.resize(n);
    dec_param.clusts_CNs.resize(n);
    dec_param.clusts_VNs.resize(n);
    dec_param.coefs_id.resize(n);

    dec_param.Roots_V[n].resize(1U << n, false);

    for (uint16_t l = 0; l < n; l++)
    {
        dec_param.Roots_V[l].resize(1U << l, false);
        dec_param.Roots_indices[l].resize(pow(2, l));
        dec_param.clusts_CNs[l].resize(pow(2, l));
        dec_param.clusts_VNs[l].resize(pow(2, l));
        dec_param.coefs_id[l].resize(pow(2, l));

        for (uint16_t s = 0; s < dec_param.Roots_V[l].size(); s++)
        {
            uint16_t sz1 = N >> (l + 1U), sz2 = sz1 << 1U;
            dec_param.clusts_CNs[l][s].resize(sz1);
            dec_param.clusts_VNs[l][s].resize(sz1);
            dec_param.coefs_id[l][s].resize(sz1);
            for (uint16_t t = 0; t < sz1; ++t)
            {
                dec_param.clusts_CNs[l][s][t] = s * sz2 + (t << 1U) - (t % sz1);
                dec_param.clusts_VNs[l][s][t] = dec_param.clusts_CNs[l][s][t] + (N >> (l + 1U));
                dec_param.coefs_id[l][s][t] = dec_param.clusts_VNs[l][s][t] - (s + 1) * sz1;
            }
        }
        for (uint16_t s = 0; s < dec_param.Roots_indices[l].size(); s++)
        {

            uint16_t sz1 = N >> l;
            dec_param.Roots_indices[l][s].resize(sz1);
            for (uint16_t t = 0; t < sz1; ++t)
                dec_param.Roots_indices[l][s][t] = s * sz1 + t;
        }
    }
    frozen_lay_pos(dec_param, dec_param.ufrozen, dec_param.clst_frozen);

    q = code_param.q;
    p = code_param.p;

    dec_param.ucap.resize(n + 1, vector<uint16_t>(N, dec_param.MxUS));
    dec_param.ucap[n].assign(N, dec_param.frozen_val);

    vector<vector<decoder_t>> L(n + 1, vector<decoder_t>(N));
    for (int i = 0; i <= n; i++)
        for (int j = 0; j < N; j++)
            L[i][j] = decoder_t(vector<softdata_t>(nm), vector<uint16_t>(nm));

    vector<vector<softdata_t>> chan_LLR(N, vector<softdata_t>(q, 0));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < q; j++)
            chan_LLR[i][j] = (softdata_t)LLR_in[i * q + j];
    }
    LLR_sort(chan_LLR, nm, L[0], dec_param.offset);
    vector<uint16_t> info_sec_rec(K, dec_param.MxUS);
    decode_SC(dec_param, table.ADDGF, table.MULGF, table.DIVGF, L, info_sec_rec);
}

#include "mex.h"
#include <vector>
#include <cstdint>

// Function declaration
void NBPolarDecoder(const std::vector<double> &LLR_in,
                    const uint16_t K,
                    const std::vector<uint16_t> &reliab_sequence,
                    const uint16_t nm,
                    const double ofst,
                    std::vector<uint16_t> &KSYMB);

// Entry point for the MEX file
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6)
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidNumInputs", "Six input arguments required.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidNumOutputs", "One output required (KSYMB).");

    // Argument 0: LLR_in (vector<double>)
    size_t len1 = mxGetNumberOfElements(prhs[0]);
    if (!mxIsDouble(prhs[0]))
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidLLR", "LLR_in must be of type double.");
    double *LLR_ptr = mxGetPr(prhs[0]);
    std::vector<double> LLR_in(LLR_ptr, LLR_ptr + len1);

    // Argument 1: K (uint16 scalar)
    if (!mxIsUint16(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidK", "K must be a scalar uint16.");
    uint16_t K = *(uint16_t *)mxGetData(prhs[1]);

    // Argument 2: reliab_sequence (vector<uint16_t>)
    if (!mxIsUint16(prhs[2]))
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidReliab", "reliab_sequence must be uint16.");
    size_t len2 = mxGetNumberOfElements(prhs[2]);
    uint16_t *reliab_ptr = (uint16_t *)mxGetData(prhs[2]);
    std::vector<uint16_t> reliab_sequence(reliab_ptr, reliab_ptr + len2);

    // Argument 3: nm (uint16 scalar)
    if (!mxIsUint16(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidNM", "nm must be a scalar uint16.");
    uint16_t nm = *(uint16_t *)mxGetData(prhs[3]);

    // Argument 4: ofst (double scalar)
    if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1)
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidOfst", "ofst must be a scalar double.");
    double ofst = mxGetScalar(prhs[4]);

    // Argument 5: initial KSYMB (vector<uint16_t> - can be empty)
    if (!mxIsUint16(prhs[5]))
        mexErrMsgIdAndTxt("NBPolarDecoder_mex:invalidKSYMB", "KSYMB must be uint16.");
    uint16_t *ksymb_ptr = (uint16_t *)mxGetData(prhs[5]);
    size_t len3 = mxGetNumberOfElements(prhs[5]);
    std::vector<uint16_t> KSYMB(ksymb_ptr, ksymb_ptr + len3);

    // Call the decoder function
    NBPolarDecoder(LLR_in, K, reliab_sequence, nm, ofst, KSYMB);

    // Return KSYMB to MATLAB
    plhs[0] = mxCreateNumericMatrix(1, KSYMB.size(), mxUINT16_CLASS, mxREAL);
    uint16_t *out_ptr = (uint16_t *)mxGetData(plhs[0]);
    std::copy(KSYMB.begin(), KSYMB.end(), out_ptr);
}
