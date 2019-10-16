//
// ij_count_diffs.c
//
// This function should conform to GNU/Intel C and Fortran interoperability standards
// and be callable from Fortran code.
//
// For two sorted (ascending order) integer vectors of dimension N, determine how many
// elements are NOT contained in both vectors.  Not that the two vectors MUST be sorted
// in ascending order for this function to work properly.
//
// The max_count limits the number of differences noted before returning, from 1 up to
// N.  If max_count <= 0, N is implied.
//
// This function is (generally speaking) complexity O(max_count) with a worst-case O(N).
//
// Build Options
// ~~~~~~~~~~~~~
// If your Fortran code is built with 64-bit integers (e.g. -i8) be sure to define the
// FORTRAN_INTEGER_64 macro.
//
//
// Vector Acceleration (AVX/AVX2/AVX512)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The function includes hand-optimized pre-check loops that utilize AVX/AVX2/AVX512
// vector instructions to (maybe?) accelerate the initial discovery of mismatched
// elements.
//
// Vector acceleration can be avoided by defining the AVOID_AVX macro.
//
// Under AVX512, the 256-bit API can be utilized by defining the AVX512_USE_256BIT
// macro.
//


#include <emmintrin.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdint.h>

#ifdef FORTRAN_INTEGER_64
    typedef int64_t             finteger;

    #ifdef AVX512_USE_256BIT
        #define AVX512_VEC_T    __m256i
        #define AVX512_LOAD_FN  _mm256_load_epi64
        #define AVX512_MASK_T   __mmask8
        #define AVX512_CMP_FN   _mm256_cmpeq_epi64_mask
    #else
        #define AVX512_VEC_T    __m512i
        #define AVX512_LOAD_FN  _mm512_load_epi64
        #define AVX512_MASK_T   __mmask8
        #define AVX512_CMP_FN   _mm512_cmpeq_epi64_mask
    #endif
    #define AVX2_LOAD_FN        _mm256_load_si256
    #define AVX_LOAD_FN         _mm_load_si128
    #define AVX2_CMP_FN         _mm256_cmpeq_epi64
    #define AVX_CMP_FN          _mm_cmpeq_epi64
#else
    typedef int32_t             finteger;

    #ifdef AVX512_USE_256BIT
        #define AVX512_VEC_T    __m256i
        #define AVX512_LOAD_FN  _mm256_load_epi32
        #define AVX512_MASK_T   __mmask8
        #define AVX512_CMP_FN   _mm256_cmpeq_epi32_mask
    #else
        #define AVX512_VEC_T    __m512i
        #define AVX512_LOAD_FN  _mm512_load_epi32
        #define AVX512_MASK_T   __mmask16
        #define AVX512_CMP_FN   _mm512_cmpeq_epi32_mask
    #endif
    #define AVX2_LOAD_FN        _mm256_load_si256
    #define AVX_LOAD_FN         _mm_load_si128
    #define AVX512_CMP_FN       _mm512_cmpeq_epi32_mask
    #define AVX2_CMP_FN         _mm256_cmpeq_epi32
    #define AVX_CMP_FN          _mm_cmpeq_epi32
#endif

//

finteger
ij_count_diffs_(
    finteger         *i,
    finteger         *j,
    finteger         *N,
    finteger         *max_count
)
{
    finteger         ii, ji, diffs;
    finteger         diff_max = *max_count;
    finteger         n = *N;

#ifdef DEBUG
    printf("Enter ij_count_diffs_(%p, %p, %d, %d)\n", i, j, n, diff_max);
#endif

#if ! defined(AVOID_AVX)
        // Implicit max_count?
        if ( diff_max <= 0 ) diff_max = n;

        // Initialize indices and difference counter:
        ii = 0;
        diffs = 0;

    #if defined(__AVX512F__) && defined(__AVX512VL__)
        //
        // AVX512
        //
        // This is the most capable ISA that makes uses the least code
        // outside the vector extensions.  It defaults to 512-bit but
        // the AVX512_USE_256BIT can override that and drop back to
        // 256-bit vectors.
        //
        // Assumes the incoming arrays are properly aligned.
        //
        finteger         rows = n / (sizeof(AVX512_VEC_T)/sizeof(finteger));

    #ifdef DEBUG
        printf("AVX512 Pre-check (rows = %d, %d-bit vectors)\n", rows, (sizeof(AVX512_VEC_T) * 8));
    #endif
        if ( rows > 0 ) {
            finteger        *iPtr = (finteger*)i;
            finteger        *jPtr = (finteger*)j;
            AVX512_VEC_T    I, J;
            AVX512_MASK_T   DIJ;
            finteger        row = 0;

            while ( row < rows ) {
                I = AVX512_LOAD_FN(iPtr);
                J = AVX512_LOAD_FN(jPtr);
                DIJ = AVX512_CMP_FN(I, J);
                if ( ! DIJ ) {
                    // Find the first non-zero bit:
                    while ( DIJ && ((DIJ & 1) == 0) ) {
                        ii++;
                        DIJ >>= 1;
                    }
                    break;
                }
                iPtr += (sizeof(AVX512_VEC_T)/sizeof(finteger)); jPtr += (sizeof(AVX512_VEC_T)/sizeof(finteger));
                ii += (sizeof(AVX512_VEC_T)/sizeof(finteger));
                row++;
            }
            if ( ii == n ) return 0;
        }
    #elif defined(__AVX2__)
        //
        // AVX2 (AVX256)
        //
        // Assumes the incoming arrays are properly aligned.  Performance
        // suffers in the check of the invidual elements of the vector
        // for the first mismatch.
        //
        finteger         rows = n / (sizeof(__m256i)/sizeof(finteger));

    #ifdef DEBUG
        printf("AVX2 Pre-check (rows = %d, %d-bit vectors)\n", rows, (sizeof(__m256i)*8));
    #endif
        if ( rows > 0 ) {
            __m256i     *iPtr = (__m256i*)i;
            __m256i     *jPtr = (__m256i*)j;
            __m256i     DIJ;
            finteger    *DIJPtr = (finteger*)&DIJ;
            finteger    row = 0;

            while ( row < rows ) {
                DIJ = AVX2_CMP_FN(*iPtr++, *jPtr++);
                if ( ! DIJPtr[0] ) break;
                ii++;
                if ( ! DIJPtr[1] ) break;
                ii++;
                if ( ! DIJPtr[2] ) break;
                ii++;
                if ( ! DIJPtr[3] ) break;
                ii++;
    #ifndef FORTRAN_INTEGER_64
                if ( ! DIJPtr[4] ) break;
                ii++;
                if ( ! DIJPtr[5] ) break;
                ii++;
                if ( ! DIJPtr[6] ) break;
                ii++;
                if ( ! DIJPtr[7] ) break;
                ii++;
    #endif
                row++;
            }
            if ( ii == n ) return 0;
        }
    #elif defined(__AVX__)
        //
        // AVX
        //
        // Assumes the incoming arrays are properly aligned.  Performance
        // suffers in the check of the invidual elements of the vector
        // for the first mismatch.
        //
        finteger         rows = n / (sizeof(__m128i)/sizeof(finteger));

    #ifdef DEBUG
        printf("AVX Pre-check (rows = %d, %d-bit vectors)\n", rows, (sizeof(__m128i)*8));
    #endif
        if ( rows > 0 ) {
            __m128i     *iPtr = (__m128i*)i;
            __m128i     *jPtr = (__m128i*)j;
            __m128i     DIJ;
            finteger    *DIJPtr = (finteger*)&DIJ;
            finteger    row = 0;

            while ( row < rows ) {
                DIJ = AVX_CMP_FN(*iPtr++, *jPtr++);
                if ( ! DIJPtr[0] ) break;
                ii++;
                if ( ! DIJPtr[1] ) break;
                ii++;
    #ifndef FORTRAN_INTEGER_64
                if ( ! DIJPtr[2] ) break;
                ii++;
                if ( ! DIJPtr[3] ) break;
                ii++;
    #endif
            }
            if ( ii == n ) return 0;
        }
    #endif
#endif

#ifdef DEBUG
    printf("ii = %d\n", ii);
#endif

    // Resume at the ii-th index (or start there if there was no vectorized
    // code pre-check:
    ji = ii;
    while ((diffs < diff_max) && (ii < n) && (ji < n)) {
        int     dij = i[ii] - j[ji];

        if ( dij == 0 ) {
            ii++;
            ji++;
        } else if ( dij < 0 ) {
            // this means i(ii) is NOT in j, so we move ahead on i
            // but NOT j and note the difference:
            ii++;
            diffs++;
        } else {
            // this means j(ji) is NOT in i, so we move ahead on j
            // but NOT i and note the difference:
            ji++;
            diffs++;
        }
    }
    // If we haven't yet hit diff_max, then we can absorb however
    // many indices were not traversed on each vector:
    if ((diffs < diff_max) && (ii < n)) {
        // Additional diffs on i are (n - ii)
        diffs = diffs + (n - ii);
        if ( diffs > diff_max ) diffs = diff_max;
    }
    if ((diffs < diff_max) && (ji < n)) {
        // Additional diffs on j are (n - ji)
        diffs = diffs + (n - ji);
        if ( diffs > diff_max ) diffs = diff_max;
    }

    return diffs;
}

//

#ifdef __BUILD_TEST__

#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>

volatile sig_atomic_t  is_running = 1;

void
alarm_raised(
    int     signum
)
{
    is_running = 0;
}

int
main()
{
    finteger            N = 56;
    finteger            i[56] = {
                            41, 42, 43, 44, 45, 46, 47, 48,
                            49, 50, 51, 52, 53, 54, 55, 56,
                            57, 58, 59, 60, 61, 62, 63, 64,
                            41, 42, 43, 44, 45, 46, 47, 48,
                            49, 50, 51, 52, 53, 54, 55, 56,
                            57, 58, 59, 60, 61, 62, 63, 64,
                            65, 66, 67, 68, 69, 70, 71, 72
                        };
    finteger            j[56] = {
                            41, 42, 43, 44, 45, 46, 47, 48,
                            49, 50, 51, 52, 53, 54, 55, 56,
                            57, 58, 59, 60, 61, 62, 63, 64,
                            41, 42, 43, 44, 45, 46, 47, 48,
                            49, 50, 51, 52, 53, 54, 55, 56,
                            57, 58, 59, 60, 61, 62, 63, 64,
                            66, 67, 68, 69, 70, 71, 72, 73
                        };
    finteger            ii, D, call_count = 0;
    struct rusage       start, end;
    double              seconds;

    signal(SIGALRM, alarm_raised);
    alarm(10);
    getrusage(RUSAGE_SELF, &start);
    while ( is_running ) {
        D = ij_count_diffs_(i, j, &N, &N);
        call_count++;
    }
    getrusage(RUSAGE_SELF, &end);

    // Calculate the number of calls per second:
    seconds = (end.ru_utime.tv_sec - start.ru_utime.tv_sec) +  (end.ru_stime.tv_sec - start.ru_stime.tv_sec);
    seconds += 1e-6 * ((end.ru_utime.tv_usec - start.ru_utime.tv_usec) +  (end.ru_stime.tv_usec - start.ru_stime.tv_usec));

    printf("%24lld calls / %-18.9f s = %24.1f calls/sec\n", call_count, seconds, (double)call_count / seconds);

    return 0;
}

#endif /* __BUILD_TEST__ */
