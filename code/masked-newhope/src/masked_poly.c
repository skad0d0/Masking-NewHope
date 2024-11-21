#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"

// convert bytes to masked poly
void poly_masked_frombytes(masked_poly *masked_r, const unsigned char *a)
{
   int i, m;
   poly *r;
   for (m = 0; m < (NEWHOPE_MASKING_ORDER + 1); m++)
   {
    r = &((masked_r->poly_shares)[m]);
    for (i = 0; i < NEWHOPE_N/4; i++)
    {
        r->coeffs[4*i+0] = (a[7*i + (m * NEWHOPE_POLYBYTES) + 0] >> 0) | 
                           (((uint16_t)a[7*i + (m * NEWHOPE_POLYBYTES) + 1] & 0x3f) << 8);

        r->coeffs[4*i+1] = (a[7*i + (m * NEWHOPE_POLYBYTES) + 1] >> 6) | 
               (((uint16_t)a[7*i + (m * NEWHOPE_POLYBYTES) + 2]) << 2) | 
               (((uint16_t)a[7*i + (m * NEWHOPE_POLYBYTES) + 3] & 0x0f) << 10);

        r->coeffs[4*i+2] = (a[7*i + (m * NEWHOPE_POLYBYTES) + 3] >> 4) | 
                (((uint16_t)a[7*i+ (m * NEWHOPE_POLYBYTES) + 4]) << 4) | 
                (((uint16_t)a[7*i+ (m * NEWHOPE_POLYBYTES) + 5] & 0x03) << 12);

        r->coeffs[4*i+3] = (a[7*i + (m * NEWHOPE_POLYBYTES) + 5] >> 2) | 
                (((uint16_t)a[7*i+ (m * NEWHOPE_POLYBYTES) + 6]) << 6);
    }
   }
}

// masked poly addition
// masked_r = masked_a + masked_b
void poly_masked_add(masked_poly *masked_r, const masked_poly *masked_a, const masked_poly *masked_b)
{
    int i,m;
    poly *r, *a, *b;
    for (m = 0; m < (NEWHOPE_MASKING_ORDER + 1); m++)
    {
        r = &((masked_r->poly_shares)[m]);
        a = &((masked_a->poly_shares)[m]);
        b = &((masked_b->poly_shares)[m]);
        for (i = 0; i < NEWHOPE_N; i++)
            r->coeffs[i] = (a->coeffs[i] + b->coeffs[i]) % NEWHOPE_Q;
    }
}

// masked poly subtraction
// masked_r = masked_a - masked_b
void poly_masked_sub(masked_poly *masked_r, const masked_poly *masked_a, const masked_poly *masked_b)
{
    int i,m;
    poly *r, *a, *b;
    for (m = 0; m < (NEWHOPE_MASKING_ORDER + 1); m++)
    {
        r = &((masked_r->poly_shares)[m]);
        a = &((masked_a->poly_shares)[m]);
        b = &((masked_b->poly_shares)[m]);
        for (i = 0; i < NEWHOPE_N; i++)
            r->coeffs[i] = (a->coeffs[i] + 3*NEWHOPE_Q - b->coeffs[i]) % NEWHOPE_Q;
    }
}

void poly_masked_mul_pointwise(masked_poly *masked_r, const masked_poly *masked_a, const masked_poly *masked_b)
{
    int i,m;
    uint16_t t;
    poly *r, *a, *b;
    for (m = 0; m < (NEWHOPE_MASKING_ORDER + 1); m++)
    {
        r = &((masked_r->poly_shares)[m]);
        a = &((masked_a->poly_shares)[m]);
        b = &((masked_b->poly_shares)[m]);
        for (i = 0; i < NEWHOPE_N; i++)
        {
            t            = montgomery_reduce(3186*b->coeffs[i]);
            r->coeffs[i] = montgomery_reduce(a->coeffs[i] * t);
        }
    }
}

void poly_halfmasked_mul_pointwise(masked_poly *masked_r, const poly *a, const masked_poly *masked_b)
{
    int i,m;
    uint16_t t;
    poly *r, *b;
    for (m = 0; m < (NEWHOPE_MASKING_ORDER + 1); m++)
    {
        r = &((masked_r->poly_shares)[m]);
        // a = &((masked_a->poly_shares)[m]);
        b = &((masked_b->poly_shares)[m]);
        for (i = 0; i < NEWHOPE_N; i++)
        {
            t            = montgomery_reduce(3186*b->coeffs[i]);
            r->coeffs[i] = montgomery_reduce(a->coeffs[i] * t);
        }
    }
}


void poly_masked_ntt(masked_poly *masked_r)
{
    int m;
    poly *r;
    for (m = 0; m < (NEWHOPE_MASKING_ORDER + 1); m++)
    {
        r = &((masked_r->poly_shares)[m]);
        mul_coefficients(r->coeffs, gammas_bitrev_montgomery);
        ntt((uint16_t *)r->coeffs, gammas_bitrev_montgomery);
    }
}

void poly_masked_invntt(masked_poly *masked_r)
{
    int m;
    poly *r;
    for (m = 0; m < (NEWHOPE_MASKING_ORDER + 1); m++)
    {
        r = &((masked_r->poly_shares)[m]);
        bitrev_vector(r->coeffs);
        ntt((uint16_t *)r->coeffs, omegas_inv_bitrev_montgomery);
        mul_coefficients(r->coeffs, gammas_inv_montgomery);
    }
}