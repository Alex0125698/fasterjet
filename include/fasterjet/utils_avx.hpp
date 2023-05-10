/**
 * \file utils_avx.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include "precompiledheader.hpp"
#include <immintrin.h>

namespace AVX256
{
   constexpr size_t DX = 8;
   using vecXf = __m256;
   using vecXd = __m256d;
   using vecXi = __m256i;
   // constexpr auto _mm256_load_ps = _mm256__mm256_load_ps;
   // constexpr auto _mm256_set1_ps = _mm256__mm256_set1_ps;
   // constexpr auto _mm256_div_ps = _mm256__mm256_div_ps;
   // constexpr auto setzero_si = _mm256_setzero_si256;
   // constexpr auto _mm256_setr_epi32 = _mm256__mm256_setr_epi32;
   // constexpr auto set_ps = _mm256_set_ps;
   // constexpr auto _mm256_set1_epi32 = _mm256__mm256_set1_epi32;
   // constexpr auto _mm256_add_epi32 = _mm256__mm256_add_epi32;
   // constexpr auto _mm256_sub_ps = _mm256__mm256_sub_ps;
   // constexpr auto _mm256_andnot_ps = _mm256_andnot_ps;
   // constexpr auto _mm256_cmpeq_epi32 = _mm256__mm256_cmpeq_epi32;
   // constexpr auto _mm256_min_ps = _mm256__mm256_min_ps;
   // constexpr auto max_ps = _mm256_max_ps;
   // constexpr auto _mm256_mul_ps = _mm256_mul_ps;
   // constexpr auto _mm256_fmadd_ps = _mm256__mm256_fmadd_ps;
   // constexpr auto cmp_ps = _mm256_cmp_ps;
   // constexpr auto blendv_epi8 = _mm256_blendv_epi8;
   // constexpr auto blendv_ps = _mm256_blendv_ps;
   // constexpr auto store_ps = _mm256_store_ps;
   // constexpr auto store_si = _mm256_store_si256;
   // constexpr auto _mm256_castsi256_ps = _mm256_castsi256_ps;
   // constexpr auto _mm256_castps_si256 = _mm256__mm256_castps_si256256;
   // constexpr auto maskstore_ps = _mm256_maskstore_ps;
   // constexpr auto maskstore_epi32 = _mm256_maskstore_epi32;
}

// namespace AVX512
// {
//    constexpr size_t DX = 16;
//    using vecXf = __m512;
//    using vecXd = __m512d;
//    using vecXi = __m512i;
//    constexpr auto _mm256_load_ps = _mm512__mm256_load_ps;
//    constexpr auto _mm256_set1_ps = _mm512__mm256_set1_ps;
//    constexpr auto _mm256_div_ps = _mm512__mm256_div_ps;
//    constexpr auto setzero_si = _mm512_setzero_si512;
//    constexpr auto _mm256_setr_epi32 = _mm512__mm256_setr_epi32;
//    constexpr auto _mm256_set1_epi32 = _mm512__mm256_set1_epi32;
//    constexpr auto _mm256_add_epi32 = _mm512__mm256_add_epi32;
//    constexpr auto _mm256_sub_ps = _mm512__mm256_sub_ps;
//    constexpr auto _mm256_andnot_ps = _mm512__mm256_andnot_ps;
//    // constexpr auto _mm256_cmpeq_epi32 = _mm512__mm256_cmpeq_epi32;
//    constexpr auto _mm256_min_ps = _mm512__mm256_min_ps;
//    constexpr auto _mm256_mul_ps = _mm512__mm256_mul_ps;
//    constexpr auto _mm256_fmadd_ps = _mm512__mm256_fmadd_ps;
//    // constexpr auto cmp_ps = _mm512_cmp_ps;
//    // constexpr auto blendv_epi8 = _mm512_blendv_epi8;
//    // constexpr auto blendv_ps = _mm512_blendv_ps;
//    constexpr auto store_ps = _mm512_store_ps;
//    // constexpr auto store_si = _mm512_store_si512;
//    constexpr auto _mm256_castsi256_ps = _mm512_castsi512_ps;
//    constexpr auto _mm256_castps_si256 = _mm512__mm256_castps_si256512;
// }
