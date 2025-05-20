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

void NBPolarEncode(    const std::vector<uint16_t>& KSYMB,
    const std::vector<uint16_t>& reliab_sequence,
    const uint16_t K,
    const uint16_t N,
    std::vector<uint16_t>& NSYMB)
{

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

    table_GF table;
    // LoadCode(code_param, EbN0);
    LoadTables(code_param, table, GF_polynom_primitive.data());

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
    CCSK_seq ccsk_seq;
    vector<vector<uint16_t>> CCSK_rotated_codes(q, vector<uint16_t>());
    if (code_param.sig_mod == "CCSK_BIN")
        create_ccsk_rotated_table(ccsk_seq.CCSK_bin_seq[code_param.p - 2], ccsk_seq.CCSK_bin_seq[code_param.p - 2].size(), CCSK_rotated_codes);
    else if (code_param.sig_mod == "CCSK_NB")
        create_ccsk_rotated_table(ccsk_seq.CCSK_GF_seq[code_param.p - 2], ccsk_seq.CCSK_GF_seq[code_param.p - 2].size(), CCSK_rotated_codes);

    vector<vector<vector<int16_t>>> hst1(n, vector<vector<int16_t>>(N, vector<int16_t>(dec_param.nm, 0)));

    q = code_param.q;
    p = code_param.p;
    vector<vector<softdata_t>> bin_mod_dict;
    if (code_param.sig_mod == "CCSK_BIN")
    {
        bin_mod_dict.resize(q, vector<softdata_t>(q, 0));

        for (int i = 0; i < q; i++)
            for (int j = 0; j < q; j++)
                bin_mod_dict[i][j] = (CCSK_rotated_codes[i][j] == 0) ? 1 : -1;
    }
    else if (code_param.sig_mod == "BPSK")
    {
        bin_mod_dict.resize(q, vector<softdata_t>(p, 0));
        for (int i = 0; i < q; i++)
            for (int j = 0; j < p; j++)
                bin_mod_dict[i][j] = (table.BINDEC[i][j] == 0) ? 1 : -1;
    }

    dec_param.ucap.resize(n + 1, vector<uint16_t>(N, dec_param.MxUS));
    dec_param.ucap[n].assign(N, dec_param.frozen_val);

    vector<vector<uint16_t>> KBIN(K);
    unsigned base_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(0 + base_seed);

    std::uniform_int_distribution<int> unif_dist(0, q - 1);
    for (int i = 0; i < K; i++)
        dec_param.ucap[dec_param.n][dec_param.reliab_sequence[i]] = KSYMB[i];
    Encoder(table.ADDGF, table.MULGF, dec_param.polar_coeff, dec_param.ucap, NSYMB);

    // vector<vector<decoder_t>> L(n + 1, vector<decoder_t>(N));
    // for (int i = 0; i <= n; i++)
    //     for (int j = 0; j < N; j++)
    //         L[i][j] = decoder_t(vector<softdata_t>(nm), vector<uint16_t>(nm));

    // vector<vector<softdata_t>> chan_LLR(N, vector<softdata_t>(q, 0));
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < q; j++)
    //         if (j == NSYMB[i])
    //             chan_LLR[i][j] = 0;
    //         else
    //             chan_LLR[i][j] = 1e5;
    // }
    // LLR_sort(chan_LLR, nm, L[0], dec_param.offset);
    // vector<uint16_t> info_sec_rec(K, dec_param.MxUS);
    // decode_SC(dec_param, table.ADDGF, table.MULGF, table.DIVGF, L, info_sec_rec);
}
