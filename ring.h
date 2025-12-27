#ifndef RING_H
#define RING_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct {
    size_t words;
    uint32_t *w;
} Ring;

void zero(Ring *a);
void copy(Ring *dst, const Ring *src);
void normalize(Ring *a);

void ring_alloc(Ring *a, size_t words);
void ring_free(Ring *a);

void add_mod(Ring *r, const Ring *a, const Ring *b);
void sub_mod(Ring *r, const Ring *a, const Ring *b);
void mul_pow2_mod(Ring *a, size_t k);
void div_pow2_mod(Ring *a, size_t k);

void add_at(Ring *dst, const Ring *src, size_t offset);


bool ring_from_hex(const char *hex, Ring *out);
bool ring_to_hex(const Ring *a, char *out, size_t out_len);
bool ring_to_bitstr(const Ring *a, char *out, size_t out_len);
bool ring_from_bitstr(const char *in, Ring *out);

void fft(Ring *a, size_t N, size_t logN);
void ifft(Ring *a, size_t N, size_t logN);
void weighted_fft(Ring *a, size_t N);
void weighted_ifft(Ring *a, size_t N);

void mul_rec_full(Ring *r, const Ring *a, const Ring *b);

void ss_mul(Ring *r, const Ring *a, const Ring *b);

#endif  
