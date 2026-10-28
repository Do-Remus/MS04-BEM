#ifndef SIMD_COMMUN_HPP_INCLUDED
#define SIMD_COMMUN_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/algebre/vecteur.hpp"

#if defined(__aarch64__) || defined(_M_ARM64)
#define SIMD_NEON
#elif defined(__AVX512F__)
#define SIMD_AVX512
#elif defined(__AVX2__)
#define SIMD_AVX2
#else
#error "No supported SIMD target."
#endif

#if (defined(SIMD_AVX2) + defined(SIMD_AVX512) + defined(SIMD_NEON)) != 1
#error "Define exactly one of SIMD_AVX2, SIMD_AVX512 or SIMD_NEON"
#endif

// Vectorized work
#ifdef SIMD_NEON
#include <arm_neon.h>
#else
#include <immintrin.h>
#endif

// -----------------------------------------------------------------------------
// Types and vector width
// -----------------------------------------------------------------------------

#ifdef SIMD_AVX2

using Simd = __m256d;
using SimdIndex = __m128i;
using SimdMask = __m256d;
#define SIMD_WIDTH 4

#elif defined(SIMD_AVX512)

using Simd = __m512d;
using SimdIndex = __m256i;
using SimdMask = __mmask8;
#define SIMD_WIDTH 8

#elif defined(SIMD_NEON)

using Simd = float64x2_t;
using SimdIndex = int32x2_t;
using SimdMask = uint64x2_t;
#define SIMD_WIDTH 2

#endif

// -----------------------------------------------------------------------------
// Basic arithmetic
// -----------------------------------------------------------------------------

inline Simd simd_zero()
{
#ifdef SIMD_AVX2
    return _mm256_setzero_pd();
#elif defined(SIMD_AVX512)
    return _mm512_setzero_pd();
#else
    return vdupq_n_f64(0.0);
#endif
}

inline Simd simd_set1(const Real x)
{
#ifdef SIMD_AVX2
    return _mm256_set1_pd(x);
#elif defined(SIMD_AVX512)
    return _mm512_set1_pd(x);
#else
    return vdupq_n_f64(x);
#endif
}

inline Simd simd_load(const Real *p)
{
#ifdef SIMD_AVX2
    return _mm256_loadu_pd(p);
#elif defined(SIMD_AVX512)
    return _mm512_loadu_pd(p);
#else
    return vld1q_f64(p);
#endif
}

inline void simd_store(Real *p, const Simd x)
{
#ifdef SIMD_AVX2
    _mm256_storeu_pd(p, x);
#elif defined(SIMD_AVX512)
    _mm512_storeu_pd(p, x);
#else
    vst1q_f64(p, x);
#endif
}

inline Simd simd_add(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_add_pd(a, b);
#elif defined(SIMD_AVX512)
    return _mm512_add_pd(a, b);
#else
    return vaddq_f64(a, b);
#endif
}

inline Simd simd_sub(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_sub_pd(a, b);
#elif defined(SIMD_AVX512)
    return _mm512_sub_pd(a, b);
#else
    return vsubq_f64(a, b);
#endif
}

inline Simd simd_mul(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_mul_pd(a, b);
#elif defined(SIMD_AVX512)
    return _mm512_mul_pd(a, b);
#else
    return vmulq_f64(a, b);
#endif
}

inline Simd simd_sqrt(const Simd x)
{
#ifdef SIMD_AVX2
    return _mm256_sqrt_pd(x);
#elif defined(SIMD_AVX512)
    return _mm512_sqrt_pd(x);
#else
    return vsqrtq_f64(x);
#endif
}

// -----------------------------------------------------------------------------
// Comparisons and blend
// -----------------------------------------------------------------------------

inline SimdMask simd_cmp_lt(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_cmp_pd(a, b, _CMP_LT_OQ);
#elif defined(SIMD_AVX512)
    return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ);
#else
    return vcltq_f64(a, b);
#endif
}

inline Simd simd_blend(const Simd a, const Simd b,
                       const SimdMask masque)
{
#ifdef SIMD_AVX2
    return _mm256_blendv_pd(a, b, masque);
#elif defined(SIMD_AVX512)
    return _mm512_mask_blend_pd(masque, a, b);
#else
    // Select b when masque is true, otherwise select a.
    return vbslq_f64(masque, b, a);
#endif
}

// -----------------------------------------------------------------------------
// Integer operations
// -----------------------------------------------------------------------------

inline SimdIndex simd_double_to_int(const Simd x)
{
#ifdef SIMD_AVX2
    return _mm256_cvttpd_epi32(x);
#elif defined(SIMD_AVX512)
    return _mm512_cvttpd_epi32(x);
#else
    // Convert two doubles to signed 32-bit integers, rounding toward zero.
    alignas(16) double values[2];
    vst1q_f64(values, x);
    return vset_lane_s32(static_cast<int32_t>(values[1]), vdup_n_s32(static_cast<int32_t>(values[0])), 1);
#endif
}

inline SimdIndex simd_index_set1(const int x)
{
#ifdef SIMD_AVX2
    return _mm_set1_epi32(x);
#elif defined(SIMD_AVX512)
    return _mm256_set1_epi32(x);
#else
    return vdup_n_s32(x);
#endif
}

inline SimdIndex simd_index_mul(const SimdIndex a, const SimdIndex b)
{
#ifdef SIMD_AVX2
    return _mm_mullo_epi32(a, b);
#elif defined(SIMD_AVX512)
    return _mm256_mullo_epi32(a, b);
#else
    return vmul_s32(a, b);
#endif
}

inline SimdIndex simd_index_add(const SimdIndex a, const SimdIndex b)
{
#ifdef SIMD_AVX2
    return _mm_add_epi32(a, b);
#elif defined(SIMD_AVX512)
    return _mm256_add_epi32(a, b);
#else
    return vadd_s32(a, b);
#endif
}

// -----------------------------------------------------------------------------
// Gather
// -----------------------------------------------------------------------------

inline Simd simd_gather(const Real *base, const SimdIndex index)
{
#ifdef SIMD_AVX2
    return _mm256_i32gather_pd(base, index, 8);
#elif defined(SIMD_AVX512)
    return _mm512_i32gather_pd(index, base, 8);
#else
    // NEON has no equivalent general-purpose double gather.
    // Gather the two elements individually.
    const Real values[2] = {
        base[vget_lane_s32(index, 0)],
        base[vget_lane_s32(index, 1)]};
    return vld1q_f64(values);
#endif
}

// -----------------------------------------------------------------------------
// Horizontal sum
// -----------------------------------------------------------------------------

inline Real simd_horizontal_sum(const Simd x)
{
#ifdef SIMD_AVX2
    const __m128d low = _mm256_castpd256_pd128(x);
    const __m128d high = _mm256_extractf128_pd(x, 1);
    const __m128d sum = _mm_add_pd(low, high);

    return _mm_cvtsd_f64(sum) + _mm_cvtsd_f64(_mm_unpackhi_pd(sum, sum));

#elif defined(SIMD_AVX512)
    return _mm512_reduce_add_pd(x);

#else
    return vgetq_lane_f64(x, 0) + vgetq_lane_f64(x, 1);
#endif
}

// -----------------------------------------------------------------------------
// Complex vector loading
// -----------------------------------------------------------------------------

inline Simd simd_load_complex_re(const Vecteur &x, const std::size_t j)
{
#ifdef SIMD_AVX2
    return _mm256_set_pd(
        x[j + 3].real(),
        x[j + 2].real(),
        x[j + 1].real(),
        x[j].real());

#elif defined(SIMD_AVX512)
    return _mm512_set_pd(
        x[j + 7].real(),
        x[j + 6].real(),
        x[j + 5].real(),
        x[j + 4].real(),
        x[j + 3].real(),
        x[j + 2].real(),
        x[j + 1].real(),
        x[j].real());

#else
    const Real values[2] = {
        x[j].real(),
        x[j + 1].real()};
    return vld1q_f64(values);
#endif
}

inline Simd simd_load_complex_im(const Vecteur &x, const std::size_t j)
{
#ifdef SIMD_AVX2
    return _mm256_set_pd(
        x[j + 3].imag(),
        x[j + 2].imag(),
        x[j + 1].imag(),
        x[j].imag());

#elif defined(SIMD_AVX512)
    return _mm512_set_pd(
        x[j + 7].imag(),
        x[j + 6].imag(),
        x[j + 5].imag(),
        x[j + 4].imag(),
        x[j + 3].imag(),
        x[j + 2].imag(),
        x[j + 1].imag(),
        x[j].imag());

#else
    const Real values[2] = {
        x[j].imag(),
        x[j + 1].imag()};
    return vld1q_f64(values);
#endif
}

inline Simd simd_index_to_real(const SimdIndex x)
{
#ifdef SIMD_AVX2
    return _mm256_cvtepi32_pd(x);
#elif defined(SIMD_AVX512)
    return _mm512_cvtepi32_pd(x);
#else
    // Conversion des deux entiers 32 bits en doubles.
    alignas(16) int32_t indices[2];
    Real values[2];

    vst1_s32(indices, x);

    values[0] = static_cast<Real>(indices[0]);
    values[1] = static_cast<Real>(indices[1]);

    return vld1q_f64(values);
#endif
}

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------

inline const Simd SIMD_HALF = simd_set1(0.5);
inline const Simd SIMD_UN = simd_set1(1.0);
inline const SimdIndex SIMD_DEUX_ID = simd_index_set1(2);
inline const SimdIndex SIMD_UN_ID = simd_index_set1(1);

#endif