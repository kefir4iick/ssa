#include "ring.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

static inline double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

void ring_random(Ring *a) {
    for (size_t i = 0; i < a->words; i++)
        a->w[i] = ((uint32_t)rand() << 16) ^ rand();
    a->w[a->words] = 0;
}

bool ring_equal(const Ring *a, const Ring *b) {
    if (a->words != b->words) return false;
    for (size_t i = 0; i <= a->words; i++)
        if (a->w[i] != b->w[i]) return false;
    return true;
}

double bench_mul(
    void (*mul)(Ring *, const Ring *, const Ring *),
    size_t words,
    int rounds
) {
    Ring a, b, r;
    ring_alloc(&a, words);
    ring_alloc(&b, words);
    ring_alloc(&r, 2 * words + 8);

    ring_random(&a);
    ring_random(&b);

    for (int i = 0; i < 3; i++)
        mul(&r, &a, &b);

    double t0 = now_sec();
    for (int i = 0; i < rounds; i++)
        mul(&r, &a, &b);
    double t1 = now_sec();

    ring_free(&a);
    ring_free(&b);
    ring_free(&r);

    return (t1 - t0) / rounds;
}

void test_correctness(size_t words, int tests) {
    Ring a, b, r1, r2;

    ring_alloc(&a,  words);
    ring_alloc(&b,  words);
    ring_alloc(&r1, 2 * words + 8);
    ring_alloc(&r2, 2 * words + 8);

    for (int t = 0; t < tests; t++) {
        ring_random(&a);
        ring_random(&b);

        mul_rec_full(&r1, &a, &b);
        ss_mul(&r2, &a, &b);

        if (!ring_equal(&r1, &r2)) {
            char sa[4096], sb[4096], s1[8192], s2[8192];

            ring_to_hex(&a,  sa, sizeof(sa));
            ring_to_hex(&b,  sb, sizeof(sb));
            ring_to_hex(&r1, s1, sizeof(s1));
            ring_to_hex(&r2, s2, sizeof(s2));

            printf("\nerror (words = %zu, test = %d)\n", words, t);
            printf("a  = %s\n", sa);
            printf("b  = %s\n", sb);
            printf("rec= %s\n", s1);
            printf("ss = %s\n", s2);

            abort();
        }
    }

    ring_free(&a);
    ring_free(&b);
    ring_free(&r1);
    ring_free(&r2);
}



int main(void) {
    srand(1);

    /*printf("| bits | words | mul_rec_full (ms) | ss_mul (ms) |\n");
    printf("|------|-------|-------------------|-------------|\n");

    for (size_t words = 16; words <= 73728; words *= 2) {
        size_t bits = words * 32;

        int rounds =
            (words <= 64)   ? 100 :
            (words <= 256)  ? 50  :
            (words <= 1024) ? 10  : 3;

        double t_rec = bench_mul(mul_rec_full, words, rounds);
        double t_ss  = bench_mul(ss_mul,       words, rounds);

        printf("| %5zu | %5zu | %17.3f | %11.3f |\n",
               bits,
               words,
               1000.0 * t_rec,
               1000.0 * t_ss);
    }*/
    
    printf("\nrun tests:\n");

    for (size_t words = 1; words <= 4096; words *= 2) {
        int tests =
            (words <= 32)   ? 200 :
            (words <= 256)  ? 100 :
            (words <= 1024) ? 30  : 10;

        printf("  words = %-5zu : %d tests\n", words, tests);
        test_correctness(words, tests);
    }

    printf("correct\n\n");

    
    
    const char *hx = "F1E2D3C4B5A697887766554433221100FFEEDDCCBBAA99887766554433221100123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0F1E2D3C4B5A697887766554433221100FFEEDDCCBBAA99887766";
                     
    const char *hy = "0FEDCBA98765432100112233445566778899AABBCCDDEEFF001122334455667789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF012345";

    Ring x, y, r1, r2;
    ring_alloc(&x, 32);   
    ring_alloc(&y, 32);
    ring_alloc(&r1, 64);  
    ring_alloc(&r2, 64);
    
    ring_from_hex(hx, &x);
    ring_from_hex(hy, &y);

    ss_mul(&r1, &x, &y);
    mul_rec_full(&r2, &x, &y);
    
    char buf1[1024];
    char buf2[1024];
    ring_to_hex(&r1, buf1, sizeof(buf1));
    ring_to_hex(&r2, buf2, sizeof(buf2));
    
    printf("mul result = %s\n", buf2);
    printf("ssa result = %s\n", buf1);

    ring_free(&x);
    ring_free(&y);
    ring_free(&r1);
    ring_free(&r2);

    return 0;
}
