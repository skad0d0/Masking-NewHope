#include "poly.h"
#include "gadgets.h"
#include "randombytes.h"
#include "cpapke.h"
#include "fips202.h"
#include "ccakem.h"
#include "masked_kem.h"
#include "utils.h"
#include "cpucycles.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define ITER 100
uint64_t start, stop;

void test_enc()
{
    unsigned char pk[NEWHOPE_CPAPKE_PUBLICKEYBYTES], sk[NEWHOPE_CPAPKE_SECRETKEYBYTES];
    unsigned char c[NEWHOPE_CPAPKE_CIPHERTEXTBYTES], masked_c[NEWHOPE_CPAPKE_CIPHERTEXTBYTES];
    unsigned char m[32], rec_m[32], masked_m[32 * (NEWHOPE_MASKING_ORDER+1)], masked_rec_m[32 * (NEWHOPE_MASKING_ORDER+1)];
    unsigned char seed[32], masked_seed[32 * (NEWHOPE_MASKING_ORDER+1)];
    masked_poly mskp;
    for (int i = 0; i < 32; i++)
    {
        m[i] = rand8();
        seed[i] = rand8();
    }
    printf("* Test IND-CPA PKE:\n");
    printf("---------------- unmasked -------------\n");
    printf(" plaintext : ");
    print_bitstring(m, 32);
    random_boolean_mask(masked_m, m, 32);
    random_boolean_mask(masked_seed, seed, 32);
    // key gen
    cpapke_keypair(pk, sk);
    poly skp;
    poly_frombytes(&skp, sk);

    for (int i = 0; i < NEWHOPE_N; i++)
        mskp.poly_shares[0].coeffs[i] = skp.coeffs[i];
    
    for (int k = 1; k < NEWHOPE_MASKING_ORDER+1; k++)
    {
        for (int i = 0; i < NEWHOPE_N; i++)
            mskp.poly_shares[k].coeffs[i] = 0;
    }

    cpapke_enc(c, m, pk, seed);
    printf(" ciphertext: ");
    print_bitstring(c, 32);
    // dec
    cpapke_dec(rec_m, c, sk);
    printf(" recovered : ");
    print_bitstring(rec_m, 32);

    // masked_enc
    printf("----------------  masked  -------------\n");
    printf(" plaintext : ");
    unmask_bitstring(masked_m, 32);
    cpapke_masked_enc(masked_c, masked_m, pk, masked_seed);
    printf(" ciphertext: ");
    print_bitstring(masked_c, 32);
    cpapke_masked_dec(masked_rec_m, masked_c, &mskp);
    printf(" recovered : ");
    unmask_bitstring(masked_rec_m, 32);
}

void test_ccakem()
{
    printf("\n* Test CCA KEM Dec:\n");
    unsigned char pk[NEWHOPE_CCAKEM_PUBLICKEYBYTES], sk[NEWHOPE_CCAKEM_SECRETKEYBYTES];
    unsigned char ct[NEWHOPE_CCAKEM_CIPHERTEXTBYTES];
    unsigned char sk_pke[NEWHOPE_CPAPKE_SECRETKEYBYTES], pkh[32];
    unsigned char masked_s[32 * (NEWHOPE_MASKING_ORDER + 1)], s[32];
    unsigned char ss[32], ss2[32], ss3[32 * (NEWHOPE_MASKING_ORDER + 1)];
    masked_poly mskp;
    poly skp;
    int i, k;

    crypto_kem_keypair(pk, sk);
    for (i = 0; i < NEWHOPE_CPAPKE_SECRETKEYBYTES; i++)
        sk_pke[i] = sk[i];
    poly_frombytes(&skp, sk_pke);

    for (int i = 0; i < NEWHOPE_N; i++)
        mskp.poly_shares[0].coeffs[i] = skp.coeffs[i];
    
    for (int k = 1; k < NEWHOPE_MASKING_ORDER+1; k++)
    {
        for (int i = 0; i < NEWHOPE_N; i++)
            mskp.poly_shares[k].coeffs[i] = 0;
    }

    for (i = 0; i < 32; i++)
        pkh[i] = sk[i + NEWHOPE_CPAPKE_SECRETKEYBYTES + NEWHOPE_CPAPKE_PUBLICKEYBYTES];

    for (i = 0; i < 32; i++)
        s[i] = sk[i + NEWHOPE_CPAPKE_SECRETKEYBYTES + NEWHOPE_CPAPKE_PUBLICKEYBYTES + 32];
    
    random_boolean_mask(masked_s, s, 32);

    crypto_kem_enc(ct, ss, pk);

    crypto_kem_dec(ss2, ct, sk);

    masked_crypto_kem_dec(ss3, ct, &mskp, pk, pkh, masked_s);
    
    printf("SharedSecret                   : "); print_bitstring(ss, 32);
    printf("Recovered SharedSecret unmasked: "); print_bitstring(ss2, 32);
    printf("Recovered SharedSecret   masked: "); unmask_bitstring(ss3, 32);
}

void test_bool_zero()
{
    unsigned char m[32], masked_m[32 * 4];
    int i, k;
    for (i = 0; i < 32; i++)
        m[i] = rand8();
    random_boolean_mask(masked_m, m, 32);

    int r;
    Masked x[32];
    for (i = 0; i < 32; i++)
    {
        for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
        {
            x[i].shares[k] = masked_m[i + k * 32];
        }
    }
    x[2].shares[2] = 3;
    r = newhope_bool_comp(x, m);
    printf("zero test = %d\n", r);
}

void timing_ccakem_dec()
{
    printf("\n* Benchmarks unmasked CCAKEM Dec: \n");
    unsigned char pk[NEWHOPE_CCAKEM_PUBLICKEYBYTES], sk[NEWHOPE_CCAKEM_SECRETKEYBYTES];
    unsigned char ct[NEWHOPE_CCAKEM_CIPHERTEXTBYTES];
    unsigned char ss[32], ss2[32];

    crypto_kem_keypair(pk, sk);
    crypto_kem_enc(ct, ss, pk);

    start = cpucycles();
    for (int i = 0; i < 10*ITER; i++)
    {
        crypto_kem_dec(ss2, ct, sk);
    }
    stop = cpucycles();
    printf("\n* Avg speed unmasked CCAKEM Dec: %.1f cycles.\n", (double)(stop-start)/(10*ITER));
}

void timing_ccakem_dec_masked()
{
    printf("\n* Benchmarks masked CCAKEM Dec: \n");
    printf("\n* MASKING_ORDER = %d\n", NEWHOPE_MASKING_ORDER);
    unsigned char pk[NEWHOPE_CCAKEM_PUBLICKEYBYTES], sk[NEWHOPE_CCAKEM_SECRETKEYBYTES];
    unsigned char ct[NEWHOPE_CCAKEM_CIPHERTEXTBYTES];
    unsigned char sk_pke[NEWHOPE_CPAPKE_SECRETKEYBYTES], pkh[32];
    unsigned char masked_s[32 * (NEWHOPE_MASKING_ORDER + 1)], s[32];
    unsigned char ss[32], ss2[32], ss3[32 * (NEWHOPE_MASKING_ORDER + 1)];
    masked_poly mskp;
    poly skp;
    int i, k;

    crypto_kem_keypair(pk, sk);
    for (i = 0; i < NEWHOPE_CPAPKE_SECRETKEYBYTES; i++)
        sk_pke[i] = sk[i];
    poly_frombytes(&skp, sk_pke);

    for (int i = 0; i < NEWHOPE_N; i++)
        mskp.poly_shares[0].coeffs[i] = skp.coeffs[i];
    
    for (int k = 1; k < NEWHOPE_MASKING_ORDER+1; k++)
    {
        for (int i = 0; i < NEWHOPE_N; i++)
            mskp.poly_shares[k].coeffs[i] = 0;
    }

    for (i = 0; i < 32; i++)
        pkh[i] = sk[i + NEWHOPE_CPAPKE_SECRETKEYBYTES + NEWHOPE_CPAPKE_PUBLICKEYBYTES];

    for (i = 0; i < 32; i++)
        s[i] = sk[i + NEWHOPE_CPAPKE_SECRETKEYBYTES + NEWHOPE_CPAPKE_PUBLICKEYBYTES + 32];
    
    random_boolean_mask(masked_s, s, 32);

    crypto_kem_enc(ct, ss, pk);

    start = cpucycles();
    for (int i = 0; i < ITER; i++)
    {
        masked_crypto_kem_dec(ss3, ct, &mskp, pk, pkh, masked_s);
    }
    stop = cpucycles();
    printf("\n* Avg speed masked CCAKEM Dec: %.1f cycles.\n", (double)(stop-start)/(ITER));
    printf("------------------------------------------------------------------------------\n");
}

void main()
{
    // test_enc();
    // test_ccakem();
    // test_bool_zero();
    timing_ccakem_dec();
    timing_ccakem_dec_masked();
}