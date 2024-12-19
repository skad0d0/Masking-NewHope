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
    unsigned char buf[128 * (NEWHOPE_MASKING_ORDER + 1)];
    int i, j, k;
    Masked a1, a2, a3, a4, b1, b2, b3, b4, y1, y2, y3, y4;
    unsigned char masked_extseed[(NEWHOPE_SYMBYTES+2) * (NEWHOPE_MASKING_ORDER+1)];
    uint64_t t;

    for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
    {
        for (i = 0; i < NEWHOPE_SYMBYTES; i++)
            masked_extseed[i + 34*k] = masked_seed[i + 32*k];
        masked_extseed[32 + 34*k] = 0;
    }
    masked_extseed[32] = nonce;


    for (j = 0; j < 8; j++) /* Generate noise in blocks of 64 coefficients */
    {
        for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
            masked_extseed[33 + 34*k] = 0; /* for each extseed[33] = j, it will be used to generate 64 coefficients*/
        masked_extseed[33] = j;

        shake256_masked(buf, 128, masked_extseed, 34);

        /* we process four coefficients at one time */
        for (i = 0; i < 64; i+=4)
        {
            for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
            {
                t = (((uint64_t) buf[2*i + 7 + k*128]) << 56) |
                    (((uint64_t) buf[2*i + 6 + k*128]) << 48) |
                    (((uint64_t) buf[2*i + 5 + k*128]) << 40) |
                    (((uint64_t) buf[2*i + 4 + k*128]) << 32) |
                    (((uint64_t) buf[2*i + 3 + k*128]) << 24) |
                    (((uint64_t) buf[2*i + 2 + k*128]) << 16) |
                    (((uint64_t) buf[2*i + 1 + k*128]) <<  8) |
                    ((uint64_t) buf[2*i + 0 + k*128]);
                
                // extract 8bits from t
                a1.shares[k] = (uint32_t) ((t >>  0) & 0xFF); b1.shares[k] = (uint32_t) ((t >>  8) & 0xFF);
                a2.shares[k] = (uint32_t) ((t >> 16) & 0xFF); b2.shares[k] = (uint32_t) ((t >> 24) & 0xFF);
                a3.shares[k] = (uint32_t) ((t >> 32) & 0xFF); b3.shares[k] = (uint32_t) ((t >> 40) & 0xFF);
                a4.shares[k] = (uint32_t) ((t >> 48) & 0xFF); b4.shares[k] = (uint32_t) ((t >> 56) & 0xFF);
            }

            // compute hw
            CBD(&a1, &b1, &y1);
            CBD(&a2, &b2, &y2);
            CBD(&a3, &b3, &y3);
            CBD(&a4, &b4, &y4);

            for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
            {
                masked_r->poly_shares[k].coeffs[j*64 + i + 0] = y1.shares[k];
                masked_r->poly_shares[k].coeffs[j*64 + i + 1] = y2.shares[k];
                masked_r->poly_shares[k].coeffs[j*64 + i + 2] = y3.shares[k];
                masked_r->poly_shares[k].coeffs[j*64 + i + 3] = y4.shares[k];
            }
        }
    }    
}

// process 2 coefficients at one time

// void poly_masked_sample(masked_poly *masked_r, const unsigned char *masked_seed, unsigned char nonce)
// {
//     unsigned char buf[128 * (NEWHOPE_MASKING_ORDER + 1)];
//     int i, j, k;
//     Masked a1, a2, a3, a4, b1, b2, b3, b4, y1, y2, y3, y4;
//     unsigned char masked_extseed[(NEWHOPE_SYMBYTES+2) * (NEWHOPE_MASKING_ORDER+1)];
//     uint32_t t;

//     for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
//     {
//         for (i = 0; i < NEWHOPE_SYMBYTES; i++)
//             masked_extseed[i + 34*k] = masked_seed[i + 32*k];
//         masked_extseed[32 + 34*k] = 0;
//     }
//     masked_extseed[32] = nonce;


//     for (j = 0; j < 8; j++) /* Generate noise in blocks of 64 coefficients */
//     {
//         for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
//             masked_extseed[33 + 34*k] = 0; /* for each extseed[33] = j, it will be used to generate 64 coefficients*/
//         masked_extseed[33] = j;

//         shake256_masked(buf, 128, masked_extseed, 34);

//         /* we process two coefficients at one time */
//         for (i = 0; i < 64; i+=2)
//         {
//             for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
//             {
//                 t = (((uint32_t) buf[2*i + 3 + k*128]) << 24) |
//                     (((uint32_t) buf[2*i + 2 + k*128]) << 16) |
//                     (((uint32_t) buf[2*i + 1 + k*128]) <<  8) |
//                     ((uint32_t)  buf[2*i + 0 + k*128]);
                
//                 // extract 8bits from t
//                 a1.shares[k] = ((t >>  0) & 0xFF); b1.shares[k] = ((t >>  8) & 0xFF);
//                 a2.shares[k] = ((t >> 16) & 0xFF); b2.shares[k] = ((t >> 24) & 0xFF);

//             }

//             // compute hw
//             CBD(&a1, &b1, &y1);
//             CBD(&a2, &b2, &y2);

//             for (k = 0; k < (NEWHOPE_MASKING_ORDER + 1); k++)
//             {
//                 masked_r->poly_shares[k].coeffs[j*64 + i + 0] = y1.shares[k];
//                 masked_r->poly_shares[k].coeffs[j*64 + i + 1] = y2.shares[k];
//             }
//         }
//     }    
// }



// msg is a boolean mask of message
void poly_masked_frommsg(masked_poly *masked_r, const unsigned char *msg)
{
    encode_message(msg, masked_r);
}



void poly_masked_tomsg(unsigned char *msg, masked_poly *masked_r)
{
    Masked ar1, ar2, bo1, bo2, t;
    int i, j, k;
    for (i = 0; i < 32*(NEWHOPE_MASKING_ORDER + 1); i++) msg[i] = 0;
    for (i = 0; i < 256; i++)
    {
        for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
        {
            ar1.shares[k] = (masked_r->poly_shares[k]).coeffs[i];
            ar2.shares[k] = (masked_r->poly_shares[k]).coeffs[i + 256];
        }
        newhope_decryption(&ar1, &bo1);
        newhope_decryption(&ar2, &bo2);

        sec_and(&bo1, &bo2, &t, 1);
        for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
        {
            msg[(i>>3) + 32*k] |= (t.shares[k]) << (i&7);
        }

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

