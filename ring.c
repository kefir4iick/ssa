#include "ring.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

void zero(Ring *a) {
    for (size_t i = 0; i <= a->words; i++)
        a->w[i] = 0;
}

void ring_alloc(Ring *a, size_t words) {
    a->words = words;
    a->w = calloc(words + 1, sizeof(uint32_t));
}

void ring_free(Ring *a) {
    if (a->w) {
        free(a->w);
        a->w = NULL;
    }
    a->words = 0;
}

void copy(Ring *dst, const Ring *src) {
    size_t min_words = dst->words < src->words ? dst->words : src->words;
    for (size_t i = 0; i < min_words; i++)
        dst->w[i] = src->w[i];
    for (size_t i = min_words; i < dst->words; i++)
        dst->w[i] = 0;
    dst->w[dst->words] = (dst->words == src->words) ? src->w[src->words] : 0;
}

static Ring ring_clone(const Ring *src) {
    Ring r;
    ring_alloc(&r, src->words);
    if (!r.w) return r;
    copy(&r, src);
    return r;
}

static inline size_t ring_bits(const Ring *a) { return a->words * 32; }

void add_mod(Ring *r, const Ring *a, const Ring *b) {
    uint64_t carry = 0;
    for (size_t i = 0; i < r->words; i++) {
        uint64_t ai = (i < a->words ? (uint64_t)a->w[i] : 0ULL);
        uint64_t bi = (i < b->words ? (uint64_t)b->w[i] : 0ULL);
        uint64_t sum = ai + bi + carry;
        r->w[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
        carry = sum >> 32;
    }
    uint64_t hi_a = (a->words == r->words ? (uint64_t)a->w[a->words] : 0ULL);
    uint64_t hi_b = (b->words == r->words ? (uint64_t)b->w[b->words] : 0ULL);
    uint64_t sum_hi = hi_a + hi_b + carry;
    r->w[r->words] = (uint32_t)sum_hi;
    carry = sum_hi >> 32;  
    if (carry) {
        uint64_t extra_carry = carry;
        for (size_t i = 0; i < r->words && extra_carry; i++) {
            uint64_t sum = (uint64_t)r->w[i] + extra_carry;
            r->w[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
            extra_carry = sum >> 32;
        }
        r->w[r->words] += extra_carry;
    }
    normalize(r);
}

void sub_mod(Ring *r, const Ring *a, const Ring *b) {
    uint64_t borrow = 0;
    for (size_t i = 0; i < r->words; i++) {
        uint64_t ai = (i < a->words ? (uint64_t)a->w[i] : 0ULL);
        uint64_t bi = (i < b->words ? (uint64_t)b->w[i] : 0ULL);
        uint64_t diff = ai - bi - borrow;
        r->w[i] = (uint32_t)(diff & 0xFFFFFFFFULL);
        borrow = (ai < bi + borrow) ? 1 : 0;
    }
    uint64_t hi_a = (a->words == r->words ? (uint64_t)a->w[a->words] : 0ULL);
    uint64_t hi_b = (b->words == r->words ? (uint64_t)b->w[b->words] : 0ULL);
    uint64_t diff_hi = hi_a - hi_b - borrow;
    r->w[r->words] = (uint32_t)(diff_hi & 0xFFFFFFFFULL);
    borrow = (hi_a < hi_b + borrow) ? 1 : 0;
    if (borrow) {
        uint64_t carry = 1;
        for (size_t i = 0; i <= r->words && carry; i++) {
            uint64_t sum = (uint64_t)r->w[i] + carry;
            r->w[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
            carry = sum >> 32;
        }
    }
    normalize(r);
}

static void left_shift_mod_2N(Ring *dst, const Ring *a, size_t r_bits) {
    size_t n = a->words;
    size_t ws = r_bits / 32;
    size_t bs = r_bits % 32;
    zero(dst);
    for (size_t i = 0; i < n; i++) {
        size_t j = i + ws;
        if (j < n) dst->w[j] = a->w[i];
        else if (j == n) dst->w[n] = a->w[i];  
    }
    if (bs) {
        uint32_t carry = 0;
        for (size_t j = 0; j <= n; j++) {  
            uint32_t cur = dst->w[j];
            uint64_t v = ((uint64_t)cur << bs) | carry;
            dst->w[j] = (uint32_t)(v & 0xFFFFFFFFULL);
            carry = (uint32_t)(v >> 32);
        }
    }
}

static void right_shift_exact(Ring *dst, const Ring *a, size_t r_bits) {
    size_t n = a->words;
    if (r_bits == 0) {
        copy(dst, a);
        return;
    }
    size_t ws = r_bits / 32;
    size_t bs = r_bits % 32;
    zero(dst);
    for (size_t i = ws; i < n; i++)
        dst->w[i - ws] = a->w[i];
    if (bs) {
        uint32_t carry = 0;
        for (size_t i = n; i-- > 0;) {
            uint32_t cur = dst->w[i];
            dst->w[i] = (cur >> bs) | carry;
            carry = cur << (32 - bs);
        }
    }
}

static void neg_mod(Ring *r, const Ring *x) {
    Ring z;
    ring_alloc(&z, x->words);
    zero(&z);
    sub_mod(r, &z, x);
    ring_free(&z);
}

void mul_pow2_mod(Ring *a, size_t k)
{
    const size_t W = a->words;
    const size_t N = W * 32;
    if (W == 0) return;

    k %= (2 * N);
    if (k == 0) return;

    bool negate = false;
    if (k >= N) {
        k -= N;
        negate = true;
    }

    const size_t ws = k >> 5;      
    const size_t bs = k & 31;      

    if (ws || bs) {
        uint32_t tmp[W + 1];
        memset(tmp, 0, sizeof(tmp));

        for (size_t i = 0; i < W; i++) {
            size_t j = i + ws;
            if (j >= W) j -= W;

            if (bs == 0) {
                tmp[j] ^= a->w[i];
            } else {
                tmp[j] ^= a->w[i] << bs;
                size_t j2 = j + 1;
                if (j2 == W) j2 = 0;
                tmp[j2] ^= a->w[i] >> (32 - bs);
            }
        }

        memcpy(a->w, tmp, W * sizeof(uint32_t));
        a->w[W] = 0;
    }

    if (negate) {
        uint64_t carry = 1;
        for (size_t i = 0; i < W; i++) {
            uint64_t v = (~a->w[i] & 0xFFFFFFFFULL) + carry;
            a->w[i] = (uint32_t)v;
            carry = v >> 32;
        }
    }
}


void div_pow2_mod(Ring *a, size_t k) {
    size_t N = ring_bits(a);
    if (N == 0) return;
    k %= (2 * N);
    if (k == 0) return;
    mul_pow2_mod(a, 2 * N - k);
}

void add_at(Ring *dst, const Ring *src, size_t offset) {
    if (offset >= dst->words) return;
    uint64_t carry = 0;
    size_t i;
    for (i = 0; i < src->words && offset + i < dst->words; i++) {
        uint64_t sum = (uint64_t)dst->w[offset + i] + src->w[i] + carry;
        dst->w[offset + i] = (uint32_t)(sum & 0xFFFFFFFFULL);
        carry = sum >> 32;
    }
    for (; carry && offset + i < dst->words; i++) {
        uint64_t sum = (uint64_t)dst->w[offset + i] + carry;
        dst->w[offset + i] = (uint32_t)(sum & 0xFFFFFFFFULL);
        carry = sum >> 32;
    }
    if (carry && offset + i == dst->words) {
        dst->w[dst->words] += carry;
    }
}

void normalize(Ring *a) {
    uint64_t hi = a->w[a->words];
    a->w[a->words] = 0;
    if (hi > 0) {
        uint64_t borrow = 0;
        Ring temp = ring_clone(a);
        for (size_t i = 0; i < a->words; i++) {
            uint64_t ai = temp.w[i];
            uint64_t bi = (i == 0 ? hi : 0);
            uint64_t diff = ai - bi - borrow;
            a->w[i] = (uint32_t)(diff & 0xFFFFFFFFULL);
            borrow = (ai < bi + borrow) ? 1 : 0;
        }
        ring_free(&temp);
        if (borrow) {
            uint64_t carry = 1;
            for (size_t i = 0; i < a->words; i++) {
                uint64_t sum = (uint64_t)a->w[i] + carry;
                a->w[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
                carry = sum >> 32;
            }
            a->w[a->words] += carry;
        }
    }
    uint32_t borrow_low = 0;
    for (size_t i = 0; i < a->words; i++) {
        uint64_t cur = a->w[i];
        uint64_t sub = borrow_low;
        uint64_t res = cur - sub;
        a->w[i] = (uint32_t)res;
        borrow_low = (res >> 63) & 1;
    }
    if (borrow_low) {
        uint64_t carry = 1;
        for (size_t i = 0; i <= a->words && carry; i++) {
            uint64_t sum = (uint64_t)a->w[i] + carry;
            a->w[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
            carry = sum >> 32;
        }
    }
}

static inline uint32_t ring_get_bit(const Ring *a, size_t bit) {
    size_t w = bit / 32, o = bit % 32;
    if (w >= a->words) return 0;
    return (a->w[w] >> o) & 1u;
}

static inline void ring_set_bit(Ring *a, size_t bit) {
    size_t w = bit / 32, o = bit % 32;
    if (w < a->words) a->w[w] |= (1u << o);
}

static int hex_val(char c) {
    if ('0' <= c && c <= '9') return c - '0';
    if ('a' <= c && c <= 'f') return c - 'a' + 10;
    if ('A' <= c && c <= 'F') return c - 'A' + 10;
    return -1;
}

bool ring_from_hex(const char *hex, Ring *out) {
    if (!hex || !out) return false;
    zero(out);
    size_t ring_bits = out->words * 32;
    size_t len = strlen(hex);
    for (size_t i = 0; i < len; i++) {
        int v = hex_val(hex[len - 1 - i]);
        if (v < 0) return false;
        for (int b = 0; b < 4; b++) {
            size_t bit = i * 4 + b;
            if (bit >= ring_bits) break;
            if (v & (1 << b)) ring_set_bit(out, bit);
        }
    }
    return true;
}

static size_t ring_effective_bits(const Ring *a) {
    if (a->w[a->words] != 0) {
        return a->words * 32 + 32 - __builtin_clz(a->w[a->words]);
    }
    for (size_t i = a->words; i-- > 0;) {
        if (a->w[i] != 0) {
            return i * 32 + 32 - __builtin_clz(a->w[i]);
        }
    }
    return 0;
}

bool ring_to_hex(const Ring *a, char *out, size_t out_len) {
    if (!a || !out) return false;
    Ring temp, neg;
    ring_alloc(&temp, a->words);
    ring_alloc(&neg, a->words);
    copy(&temp, a);
    neg_mod(&neg, a);
    size_t bits_a = ring_effective_bits(&temp);
    size_t bits_neg = ring_effective_bits(&neg);
    const Ring *chosen = (bits_neg < bits_a) ? &neg : &temp;
    size_t ring_bits = chosen->words * 32;
    size_t hex_len = (ring_bits + 3) / 4;
    if (out_len < hex_len + 1) {
        ring_free(&temp);
        ring_free(&neg);
        return false;
    }
    static const char H[] = "0123456789ABCDEF";
    memset(out, '0', hex_len);
    out[hex_len] = '\0';
    for (size_t i = 0; i < hex_len; i++) {
        uint32_t v = 0;
        for (int b = 0; b < 4; b++) {
            v |= ring_get_bit(chosen, i * 4 + b) << b;
        }
        out[hex_len - 1 - i] = H[v];
    }
    size_t s = 0;
    while (s + 1 < hex_len && out[s] == '0') s++;
    if (s) memmove(out, out + s, hex_len - s + 1);
    ring_free(&temp);
    ring_free(&neg);
    return true;
}

bool ring_to_bitstr(const Ring *a, char *out, size_t out_len) {
    (void)a; (void)out; (void)out_len;
    return false;
}

bool ring_from_bitstr(const char *in, Ring *out) {
    (void)in; (void)out;
    return false;
}

static size_t bitrev(size_t x, size_t log_n) {
    size_t r = 0;
    for (size_t i = 0; i < log_n; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

void fft(Ring *a, size_t N, size_t logN) {
    for (size_t i = 0; i < N; i++) {
        size_t j = bitrev(i, logN);
        if (j > i) {
            Ring tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }
    size_t n_bits = a[0].words * 32;
    for (size_t s = 1; s <= logN; s++) {
        size_t m = 1u << s;
        size_t half = m >> 1;
        size_t step = (2 * n_bits) / m;
        for (size_t k = 0; k < N; k += m) {
            for (size_t j = 0; j < half; j++) {
                Ring u = ring_clone(&a[k + j]);
                Ring v = ring_clone(&a[k + j + half]);
                mul_pow2_mod(&v, j * step);
                Ring x = ring_clone(&u);
                add_mod(&x, &x, &v);
                Ring y = ring_clone(&u);
                sub_mod(&y, &y, &v);
                copy(&a[k + j], &x);
                copy(&a[k + j + half], &y);
                ring_free(&u);
                ring_free(&v);
                ring_free(&x);
                ring_free(&y);
            }
        }
    }
}

void ifft(Ring *a, size_t N, size_t logN) {
    size_t n_bits = a[0].words * 32;
    for (size_t s = logN; s-- > 0;) {
        size_t m = 1u << (s + 1);
        size_t half = m >> 1;
        size_t step = (2 * n_bits) / m;
        for (size_t k = 0; k < N; k += m) {
            for (size_t j = 0; j < half; j++) {
                Ring x = ring_clone(&a[k + j]);
                Ring y = ring_clone(&a[k + j + half]);
                Ring u = ring_clone(&x);
                add_mod(&u, &u, &y);
                div_pow2_mod(&u, 1);
                Ring v = ring_clone(&x);
                sub_mod(&v, &v, &y);
                div_pow2_mod(&v, 1);
                size_t inv = (2 * n_bits - (j * step) % (2 * n_bits)) % (2 * n_bits);
                mul_pow2_mod(&v, inv);
                copy(&a[k + j], &u);
                copy(&a[k + j + half], &v);
                ring_free(&x);
                ring_free(&y);
                ring_free(&u);
                ring_free(&v);
            }
        }
    }
    for (size_t i = 0; i < N; i++) {
        size_t j = bitrev(i, logN);
        if (j > i) {
            Ring tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }
}

#define SS_BASE_WORDS 16

void mul_rec_full(Ring *r, const Ring *a, const Ring *b) {
    zero(r);
    if (a->words <= SS_BASE_WORDS || b->words <= SS_BASE_WORDS) {
        for (size_t i = 0; i < a->words; i++) {
            uint64_t carry = 0;
            for (size_t j = 0; j < b->words && i + j < r->words; j++) {
                uint64_t v = (uint64_t)a->w[i] * b->w[j] + r->w[i + j] + carry;
                r->w[i + j] = (uint32_t)v;
                carry = v >> 32;
            }
            if (i + b->words < r->words)
                r->w[i + b->words] += carry;
        }
        return;
    }
    size_t n = a->words > b->words ? a->words : b->words;
    size_t h = (n + 1) / 2;
    Ring a0, a1, b0, b1, p0, p1, p2, s0, s1;
    ring_alloc(&a0, h);
    ring_alloc(&a1, a->words > h ? a->words - h : 0);
    ring_alloc(&b0, h);
    ring_alloc(&b1, b->words > h ? b->words - h : 0);
    for (size_t i = 0; i < a->words; i++) {
        if (i < h)
            a0.w[i] = a->w[i];
        else
            a1.w[i - h] = a->w[i];
    }
    for (size_t i = 0; i < b->words; i++) {
        if (i < h)
            b0.w[i] = b->w[i];
        else
            b1.w[i - h] = b->w[i];
    }
    ring_alloc(&p0, 2 * h + 1);
    ring_alloc(&p1, a1.words + b1.words + 1);
    ring_alloc(&p2, 2 * h + 1);
    ring_alloc(&s0, h + 1);
    ring_alloc(&s1, h + 1);
    copy(&s0, &a0);
    add_mod(&s0, &s0, &a1);
    copy(&s1, &b0);
    add_mod(&s1, &s1, &b1);
    mul_rec_full(&p0, &a0, &b0);
    mul_rec_full(&p1, &a1, &b1);
    mul_rec_full(&p2, &s0, &s1);
    sub_mod(&p2, &p2, &p0);
    sub_mod(&p2, &p2, &p1);
    zero(r);
    add_at(r, &p0, 0);
    add_at(r, &p2, h);
    add_at(r, &p1, 2 * h);
    ring_free(&a0);
    ring_free(&a1);
    ring_free(&b0);
    ring_free(&b1);
    ring_free(&p0);
    ring_free(&p1);
    ring_free(&p2);
    ring_free(&s0);
    ring_free(&s1);
}

void ring_mul_mod(Ring *r, const Ring *a, const Ring *b) {
    Ring tmp;
    ring_alloc(&tmp, a->words + b->words + 1);
    mul_rec_full(&tmp, a, b);
    normalize(&tmp);
    copy(r, &tmp);
    ring_free(&tmp);
}

void weighted_fft(Ring *a, size_t N) {
    size_t n_bits = a[0].words * 32;
    for (size_t i = 0; i < N; i++) {
        size_t shift = (2 * n_bits * i) / N;
        mul_pow2_mod(&a[i], shift);
    }
    fft(a, N, __builtin_ctz(N));
}

void weighted_ifft(Ring *a, size_t N) {
    ifft(a, N, __builtin_ctz(N));
    size_t n_bits = a[0].words * 32;
    for (size_t i = 0; i < N; i++) {
        size_t shift = (2 * n_bits * i) / N;
        div_pow2_mod(&a[i], shift);
    }
}


void ss_mul(Ring *r, const Ring *a, const Ring *b) {
    size_t a_bits = a->words * 32;
    size_t b_bits = b->words * 32;
    size_t total_bits = a_bits + b_bits;
    size_t K = 1;
    size_t M_words = 0;
    while (1) {
        size_t bits_per_block = (total_bits + K - 1) / K;
        M_words = (bits_per_block + 31) / 32;
        if (K >= 1 && (M_words * 32 * K) >= total_bits) break;
        K <<= 1;
    }
    size_t M_bits = M_words * 32;
    size_t logK = __builtin_ctz(K);
    size_t inner_words = 2 * M_words + 4;
    Ring *A = calloc(K, sizeof(Ring));
    Ring *B = calloc(K, sizeof(Ring));
    for (size_t i = 0; i < K; i++) {
        ring_alloc(&A[i], inner_words);
        ring_alloc(&B[i], inner_words);
    }
    for (size_t i = 0; i < K; i++) {
        size_t offset = i * M_words;
        for (size_t w = 0; w < M_words; w++) {
            if (offset + w < a->words) A[i].w[w] = a->w[offset + w];
            if (offset + w < b->words) B[i].w[w] = b->w[offset + w];
        }
    }
    weighted_fft(A, K);
    weighted_fft(B, K);
    for (size_t i = 0; i < K; i++) {
        Ring t;
        ring_alloc(&t, 2 * inner_words + 2);
        mul_rec_full(&t, &A[i], &B[i]);
        normalize(&t);
        copy(&A[i], &t);
        ring_free(&t);
    }
    weighted_ifft(A, K);
    ring_alloc(r, a->words + b->words + 8);
    zero(r);
    for (size_t i = 0; i < K; i++) {
        add_at(r, &A[i], i * M_words);
    }
    normalize(r);
    for (size_t i = 0; i < K; i++) {
        ring_free(&A[i]);
        ring_free(&B[i]);
    }
    free(A);
    free(B);
}
