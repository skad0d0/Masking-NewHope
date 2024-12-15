#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "utils.h"
#include "ntt.h"
#include "reduce.h"
#include "fips202.h"
#include "gadgets.h"

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

// r = a - b
void poly_halfmasked_sub(masked_poly *masked_r, const masked_poly *masked_a, const poly *b)
{
    unsigned int i, m;
    poly *r;
    const poly *a;

    r = &((masked_r->poly_shares)[0]); 
    a = &((masked_a->poly_shares)[0]);
    for (i = 0; i < NEWHOPE_N; i++)
        r->coeffs[i] = (a->coeffs[i] + 3 * NEWHOPE_Q - b->coeffs[i]) % NEWHOPE_Q;

    for (m = 1; m < (NEWHOPE_MASKING_ORDER + 1); m++)
    {
        r = &((masked_r->poly_shares)[m]); 
        a = &((masked_a->poly_shares)[m]);
        for (i = 0; i < NEWHOPE_N; i++)
            r->coeffs[i] = (a->coeffs[i] + 3 * NEWHOPE_Q) % NEWHOPE_Q;
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

void poly_masked_sample(masked_poly *masked_r, const unsigned char *masked_seed, unsigned char nonce)
{
#if NEWHOPE_K != 8
#error "poly_sample in poly.c only supports k=8"
#endif
    unsigned char buf[128], a, b;
    int i, j, m;
    poly *r;
    unsigned char extseed[(NEWHOPE_SYMBYTES+2) * (NEWHOPE_MASKING_ORDER+1)];

    for (m = 0; m < (NEWHOPE_MASKING_ORDER+1); m++)
    {
        for (i = 0; i < NEWHOPE_SYMBYTES; i++)
            extseed[m*34 + i] = masked_seed[m*34 + i];
        extseed[m*34 + 32] = nonce;
        r = &((masked_r->poly_shares)[m]);
        for (i = 0; i < NEWHOPE_N/64; i++)
        {
            extseed[m*34 + 33] = i;
            shake256_masked(buf, 128, &extseed[m*34], 34);
            for (j = 0; j < 64; j++)
            {
                a = buf[2*j];
                b = buf[2*j+1];
                r->coeffs[64*i+j] = hw(a) + NEWHOPE_Q - hw(b); // to be fixed sum of binomial
            }
        } 
    }
}

// msg is a boolean mask of message
void poly_masked_frommsg(masked_poly *masked_r, const unsigned char *msg)
{
    encode_message(msg, masked_r);
}


/*
void poly_masked_tomsg(unsigned char *msg, masked_poly *masked_r)
{
    Masked ar1, ar2, bo1, bo2;
    int i, j, k;
    for (i = 0; i < 32*(NEWHOPE_MASKING_ORDER + 1); i++) msg[i] = 0;
    for (i = 0; i < 256; i++)
    {
        for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
        {
            ar1.shares[k] = (masked_r->poly_shares[k]).coeffs[i];
            ar2.shares[k] = (masked_r->poly_shares[k]).coeffs[i + 256];
            // printf("ar1.shares: %x\n", ar1.shares[k]);
            // printf("ar2.shares: %x\n", ar2.shares[k]);
            // printf("------------------------------\n");
        }
        newhope_decryption(&ar1, &bo1);
        newhope_decryption(&ar2, &bo2);
        for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
        {
            // printf("bo1.shares: %x\n", bo1.shares[k]&1);
            // printf("bo2.shares: %x\n", bo2.shares[k]&1);
            // printf("------------------------------\n");
            msg[(i>>3) + 32*k] |= (bo2.shares[k]&1) << (i&7);
        }

    }

}
*/

void poly_masked_tomsg(unsigned char *msg, masked_poly *masked_r)
{
    Masked ar1, ar2, bo1, bo2;
    unsigned char rec_m[32];
    int i, j, k, t1, t2;

    for (i = 0; i < 32; i++) rec_m[i] = 0;

    for (i = 0; i < 256; i++)
    {
        for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
        {
            ar1.shares[k] = (masked_r->poly_shares[k]).coeffs[i];
            ar2.shares[k] = (masked_r->poly_shares[k]).coeffs[i + 256];
        }

        newhope_decryption(&ar1, &bo1);
        newhope_decryption(&ar2, &bo2);
        t1 = bo1.shares[0]&1;
        t2 = bo2.shares[0]&1;

        for (k = 1; k < NEWHOPE_MASKING_ORDER + 1; k++)
        {
            t1 ^= bo1.shares[k]&1;
            t2 ^= bo2.shares[k]&1;
        }
        t1 &= t2;
        rec_m[(i>>3)] |= (t1) << (i&7);
        random_boolean_mask(msg, &rec_m);
    }

}



void unmask_poly(masked_poly* mp, poly* p){
  int16_t temp;
  for(int i=0; i < NEWHOPE_N; ++i){
    temp = 0;
    for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j){
      temp = (temp + (mp->poly_shares[j].coeffs[i]))%NEWHOPE_Q;
    }
    p->coeffs[i]=temp;
  }
}

