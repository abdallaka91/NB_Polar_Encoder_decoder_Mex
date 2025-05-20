#ifndef INIT_H
#define INIT_H

#include "struct.h"
#include <vector>
namespace PoAwN
{
    namespace init
    {
        void LoadCode(PoAwN::structures::base_code_t &code,
                      float SNR);
        void LoadTables(PoAwN::structures::base_code_t &code,
                        PoAwN::structures::table_GF &table,
                        const uint16_t *GF_polynom_primitive);
        void LoadBubblesIndcatorlists(PoAwN::structures::decoder_parameters &dec,
                                      const float SNR);
    } // namespace init
} // namespace PoAwN
#endif