#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "utils.h"
#include "ntt.h"
#include "reduce.h"
#include "fips202.h"
#include "gadgets.h"

/**
 * @brief Converts a message to a masked polynomial.
 *
 * This function takes a message `msg` and converts it into a masked polynomial `masked_r`.
 * The conversion process involves encoding the message into the polynomial shares of the masked polynomial.
 *
 * @param[out] masked_r Pointer to the masked polynomial structure where the result will be stored.
 * @param[in] msg Pointer to the input message buffer.
 */
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


/**
 * @brief Adds two masked polynomials and stores the result in a third masked polynomial.
 *
 * This function performs element-wise addition of two masked polynomials, `masked_a` and `masked_b`,
 * and stores the result in `masked_r`. The addition is performed modulo `NEWHOPE_Q`.
 *
 * @param[out] masked_r The resulting masked polynomial after addition.
 * @param[in] masked_a The first masked polynomial to be added.
 * @param[in] masked_b The second masked polynomial to be added.
 */
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


/**
 * @brief Subtracts two masked polynomials and stores the result in a third masked polynomial.
 *
 * This function performs element-wise subtraction of two masked polynomials, `masked_a` and `masked_b`,
 * and stores the result in `masked_r`. The subtraction is performed for each share of the masked polynomials.
 *
 * @param[out] masked_r The masked polynomial to store the result of the subtraction.
 * @param[in] masked_a The first masked polynomial operand.
 * @param[in] masked_b The second masked polynomial operand.
 */
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

/**
 * @brief Subtracts a polynomial from a half-masked polynomial.
 *
 * This function performs the subtraction of a polynomial `b` from a half-masked polynomial `masked_a`
 * and stores the result in `masked_r`. The subtraction is performed modulo `NEWHOPE_Q`.
 *
 * @param[out] masked_r The result of the subtraction, stored as a masked polynomial.
 * @param[in] masked_a The half-masked polynomial from which `b` is subtracted.
 * @param[in] b The polynomial to be subtracted from `masked_a`.
 */
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

/**
 * @brief Multiplies two masked polynomials pointwise and stores the result in a masked polynomial.
 *
 * This function performs a pointwise multiplication of two masked polynomials `masked_a` and `masked_b`,
 * and stores the result in the masked polynomial `masked_r`. The multiplication is performed for each
 * share of the masked polynomials.
 *
 * @param[out] masked_r Pointer to the masked polynomial where the result will be stored.
 * @param[in] masked_a Pointer to the first masked polynomial operand.
 * @param[in] masked_b Pointer to the second masked polynomial operand.
 */
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

/**
 * @brief Multiplies a polynomial with a masked polynomial pointwise and stores the result in a masked polynomial.
 *
 * This function performs a pointwise multiplication of a polynomial `a` with a masked polynomial `masked_b` and stores
 * the result in the masked polynomial `masked_r`. The multiplication is performed for each share of the masked polynomials.
 *
 * @param[out] masked_r The resulting masked polynomial after pointwise multiplication.
 * @param[in] a The input polynomial to be multiplied.
 * @param[in] masked_b The input masked polynomial to be multiplied.
 */
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


/**
 * @brief Applies the Number Theoretic Transform (NTT) to each polynomial share in a masked polynomial.
 *
 * This function iterates over each polynomial share in the masked polynomial and performs the following steps:
 * 1. Multiplies the coefficients of the polynomial by precomputed gamma values in Montgomery form.
 * 2. Applies the NTT to the polynomial coefficients.
 *
 * @param masked_r Pointer to the masked polynomial structure containing the polynomial shares.
 */
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

/**
 * @brief Perform the inverse Number Theoretic Transform (NTT) on a masked polynomial.
 *
 * This function takes a masked polynomial and applies the inverse NTT to each of its shares.
 * The process involves bit-reversing the coefficients, performing the inverse NTT using
 * precomputed omegas, and multiplying the coefficients by precomputed gamma values.
 *
 * @param masked_r Pointer to the masked polynomial to be transformed.
 */
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

/**
 * @brief Implementation of masked polynomial sampling for NewHope cryptographic scheme.
 *
 * This function generates a masked polynomial sample using a masked seed and a nonce.
 * It processes the coefficients in blocks of 64, and each block is processed four coefficients at a time.
 *
 * @param[out] masked_r       Pointer to the masked polynomial structure to store the result.
 * @param[in]  masked_seed    Pointer to the masked seed used for sampling.
 * @param[in]  nonce          Nonce value used for sampling.
 */
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




/**
 * @brief Converts a message into a masked polynomial representation.
 *
 * This function takes a message in the form of an unsigned char array and
 * encodes it into a masked polynomial structure.
 *
 * @param[out] masked_r Pointer to the masked polynomial structure where the encoded message will be stored.
 * @param[in] msg Pointer to the message to be encoded.
 */
void poly_masked_frommsg(masked_poly *masked_r, const unsigned char *msg)
{
    encode_message(msg, masked_r);
}



/**
 * @brief Converts a masked polynomial to a message.
 *
 * This function takes a masked polynomial and converts it into a message
 * by performing decryption and secure AND operations on the polynomial shares.
 *
 * @param msg Pointer to the output message buffer. The buffer should be large enough to hold the resulting message.
 * @param masked_r Pointer to the input masked polynomial structure.
 *
 * The function initializes the message buffer to zero, then iterates over each coefficient
 * of the polynomial shares, decrypts them, performs a secure AND operation, and updates the message buffer.
 */
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

/**
 * @brief Unmasks a masked polynomial.
 *
 * This function takes a masked polynomial `mp` and combines its shares to 
 * produce the original polynomial `p`. The unmasking process involves 
 * summing the coefficients of the shares and taking the result modulo `NEWHOPE_Q`.
 *
 * @param mp Pointer to the masked polynomial structure.
 * @param p Pointer to the polynomial structure where the result will be stored.
 */
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

