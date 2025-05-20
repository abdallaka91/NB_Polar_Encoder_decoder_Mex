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

void NBPolarEncode(const std::vector<uint16_t> &KSYMB,
                   const std::vector<uint16_t> &reliab_sequence,
                   std::vector<uint16_t> &NSYMB)
{
        
    uint16_t K = KSYMB.size();
    uint16_t N = reliab_sequence.size();

    NSYMB.resize(N, 0);
    string sig_mod = "CCSK_NB";
    uint16_t frozen_val = 0;
    uint16_t n = log2(N);
    uint16_t q = 64;
    uint16_t p = log2(q);
    float EbN0 = -7.5;
    softdata_t offset = 1;
    uint16_t nm = 18;
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
   

    q = code_param.q;
    p = code_param.p;
   

    dec_param.ucap.resize(n + 1, vector<uint16_t>(N, dec_param.MxUS));
    dec_param.ucap[n].assign(N, dec_param.frozen_val);

    vector<vector<uint16_t>> KBIN(K);
    for (int i = 0; i < K; i++)
        dec_param.ucap[dec_param.n][dec_param.reliab_sequence[i]] = KSYMB[i];
    Encoder(table.ADDGF, table.MULGF, dec_param.polar_coeff, dec_param.ucap, NSYMB);
}
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgIdAndTxt("NBPolarEncode_mex:input", "Two inputs required: KSYMB and reliab_sequence.");

    // Extract KSYMB
    if (!mxIsUint16(prhs[0]) || mxGetNumberOfElements(prhs[0]) < 1)
        mexErrMsgIdAndTxt("NBPolarEncode_mex:input", "KSYMB must be a non-empty uint16 vector.");
    uint16_t *KSYMB_ptr = (uint16_t *)mxGetData(prhs[0]);
    size_t K_len = mxGetNumberOfElements(prhs[0]);
    std::vector<uint16_t> KSYMB(KSYMB_ptr, KSYMB_ptr + K_len);

    // Extract reliab_sequence
    if (!mxIsUint16(prhs[1]) || mxGetNumberOfElements(prhs[1]) < 1)
        mexErrMsgIdAndTxt("NBPolarEncode_mex:input", "reliab_sequence must be a non-empty uint16 vector.");
    uint16_t *rel_ptr = (uint16_t *)mxGetData(prhs[1]);
    size_t rel_len = mxGetNumberOfElements(prhs[1]);
    std::vector<uint16_t> reliab_sequence(rel_ptr, rel_ptr + rel_len);

    // Prepare output vector
    std::vector<uint16_t> NSYMB;
    NSYMB.reserve(rel_len);

    // Call encode function
    NBPolarEncode(KSYMB, reliab_sequence, NSYMB);

    // Create MATLAB output (1 x N uint16)
    plhs[0] = mxCreateNumericMatrix(1, rel_len, mxUINT16_CLASS, mxREAL);
    uint16_t *out_ptr = (uint16_t *)mxGetData(plhs[0]);
    for (size_t i = 0; i < rel_len; ++i)
    {
        out_ptr[i] = NSYMB[i];
    }
}
