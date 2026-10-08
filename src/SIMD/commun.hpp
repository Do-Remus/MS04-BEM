#ifndef SIMD_COMMUN_HPP_INCLUDED
#define SIMD_COMMUN_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/algebre/vecteur.hpp"

// #define SIMD_AVX2
#define SIMD_AVX512

#ifdef SIMD_AVX2
using Simd = __m256d;
using SimdIndex = __m128i;
using SimdMask = __m256d;
#define SIMD_WIDTH 4
#else
using Simd = __m512d;
using SimdIndex = __m256i;
using SimdMask = __mmask8;
#define SIMD_WIDTH 8
#endif

inline Simd simd_zero()
{
#ifdef SIMD_AVX2
    return _mm256_setzero_pd();
#else
    return _mm512_setzero_pd();
#endif
}

inline Simd simd_set1(const Real x)
{
#ifdef SIMD_AVX2
    return _mm256_set1_pd(x);
#else
    return _mm512_set1_pd(x);
#endif
}

inline Simd simd_load(const Real *p)
{
#ifdef SIMD_AVX2
    return _mm256_loadu_pd(p);
#else
    return _mm512_loadu_pd(p);
#endif
}

inline void simd_store(Real *p, const Simd x)
{
#ifdef SIMD_AVX2
    _mm256_storeu_pd(p, x);
#else
    _mm512_storeu_pd(p, x);
#endif
}

inline Simd simd_add(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_add_pd(a, b);
#else
    return _mm512_add_pd(a, b);
#endif
}

inline Simd simd_sub(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_sub_pd(a, b);
#else
    return _mm512_sub_pd(a, b);
#endif
}

inline Simd simd_mul(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_mul_pd(a, b);
#else
    return _mm512_mul_pd(a, b);
#endif
}

inline Simd simd_sqrt(const Simd x)
{
#ifdef SIMD_AVX2
    return _mm256_sqrt_pd(x);
#else
    return _mm512_sqrt_pd(x);
#endif
}

inline SimdMask simd_cmp_lt(const Simd a, const Simd b)
{
#ifdef SIMD_AVX2
    return _mm256_cmp_pd(a, b, _CMP_LT_OQ);
#else
    return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ);
#endif
}

inline Simd simd_blend(const Simd a, const Simd b, const SimdMask masque)
{
#ifdef SIMD_AVX2
    return _mm256_blendv_pd(a, b, masque);
#else
    return _mm512_mask_blend_pd(masque, a, b);
#endif
}

inline SimdIndex simd_double_to_int(const Simd x)
{
#ifdef SIMD_AVX2
    return _mm256_cvttpd_epi32(x);
#else
    return _mm512_cvttpd_epi32(x);
#endif
}

inline SimdIndex simd_index_set1(const int x)
{
#ifdef SIMD_AVX2
    return _mm_set1_epi32(x);
#else
    return _mm256_set1_epi32(x);
#endif
}

inline SimdIndex simd_index_mul(const SimdIndex a, const SimdIndex b)
{
#ifdef SIMD_AVX2
    return _mm_mullo_epi32(a, b);
#else
    return _mm256_mullo_epi32(a, b);
#endif
}

inline SimdIndex simd_index_add(const SimdIndex a, const SimdIndex b)
{
#ifdef SIMD_AVX2
    return _mm_add_epi32(a, b);
#else
    return _mm256_add_epi32(a, b);
#endif
}

inline Simd simd_gather(const Real *base, const SimdIndex index)
{
#ifdef SIMD_AVX2
    return _mm256_i32gather_pd(base, index, 8);
#else
    return _mm512_i32gather_pd(index, base, 8);
#endif
}

inline Real simd_horizontal_sum(const Simd x)
{
#ifdef SIMD_AVX2
    const __m128d low = _mm256_castpd256_pd128(x);
    const __m128d high = _mm256_extractf128_pd(x, 1);
    const __m128d sum = _mm_add_pd(low, high);

    return _mm_cvtsd_f64(sum) + _mm_cvtsd_f64(_mm_unpackhi_pd(sum, sum));
#else
    return _mm512_reduce_add_pd(x);
#endif
}

inline Simd simd_load_complex_re(const Vecteur &x, const std::size_t j)
{
#ifdef SIMD_AVX2
    return _mm256_set_pd(
        x[j + 3].real(),
        x[j + 2].real(),
        x[j + 1].real(),
        x[j].real());
#else
    return _mm512_set_pd(
        x[j + 7].real(),
        x[j + 6].real(),
        x[j + 5].real(),
        x[j + 4].real(),
        x[j + 3].real(),
        x[j + 2].real(),
        x[j + 1].real(),
        x[j].real());
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
#else
    return _mm512_set_pd(
        x[j + 7].imag(),
        x[j + 6].imag(),
        x[j + 5].imag(),
        x[j + 4].imag(),
        x[j + 3].imag(),
        x[j + 2].imag(),
        x[j + 1].imag(),
        x[j].imag());
#endif
}

inline const Simd SIMD_HALF = simd_set1(0.5);
inline const SimdIndex SIMD_DEUX = simd_index_set1(2);
inline const SimdIndex SIMD_UN = simd_index_set1(1);

#endif