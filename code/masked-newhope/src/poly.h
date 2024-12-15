#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/* 
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1] 
 */
typedef struct {
  uint16_t coeffs[NEWHOPE_N];
} poly __attribute__ ((aligned (32)));

// masked poly
typedef struct 
{
  poly poly_shares[NEWHOPE_MASKING_ORDER + 1];
} masked_poly;

unsigned char hw(unsigned char a);
uint16_t coeff_freeze(uint16_t x);
uint16_t flipabs(uint16_t x);
// unmasked polynomial operation
void poly_uniform(poly *a, const unsigned char *seed);
void poly_sample(poly *r, const unsigned char *seed, unsigned char nonce);
void poly_add(poly *r, const poly *a, const poly *b);

void poly_ntt(poly *r);
void poly_invntt(poly *r);
void poly_mul_pointwise(poly *r, const poly *a, const poly *b);

void poly_frombytes(poly *r, const unsigned char *a);
void poly_tobytes(unsigned char *r, const poly *p);
void poly_compress(unsigned char *r, const poly *p);
void poly_decompress(poly *r, const unsigned char *a);

void poly_frommsg(poly *r, const unsigned char *msg);
void poly_tomsg(unsigned char *msg, const poly *x);
void poly_sub(poly *r, const poly *a, const poly *b);


// masked polynomial operation
void poly_masked_frombytes(masked_poly *masked_r, const unsigned char *a);
void poly_masked_frommsg(masked_poly *masked_r, const unsigned char *msg);
void poly_masked_tomsg(unsigned char *msg, masked_poly *masked_r); 

void poly_masked_ntt(masked_poly *masked_r);
void poly_masked_invntt(masked_poly *masked_r);

void poly_masked_mul_pointwise(masked_poly *masked_r, const masked_poly *masked_a, const masked_poly *masked_b);
void poly_halfmasked_mul_pointwise(masked_poly *masked_r, const poly *a, const masked_poly *masked_b);
void poly_masked_add(masked_poly *masked_r, const masked_poly *masked_a, const masked_poly *masked_b);
void poly_masked_sub(masked_poly *masked_r, const masked_poly *masked_a, const masked_poly *masked_b);
void poly_halfmasked_sub(masked_poly *masked_r, const masked_poly *masked_a, const poly *b);

void poly_masked_sample(masked_poly *masked_r, const unsigned char *masked_seed, unsigned char nonce);

void unmask_poly(masked_poly* mp, poly* p);
#endif
