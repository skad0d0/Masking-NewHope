#include "poly.h"
#include "gadgets.h"
#include "randombytes.h"
#include "cpapke.h"
#include "fips202.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

// test message encoding and decoding
void test_msg()
{
    unsigned char m[32*4];
    int i,k;
    for (i = 0; i < 32*4; i++)
    {
        m[i] = rand() & 0xff;
        // for (k = 0; k < 4; k++)
        //     m[i + 32*k] = i+k;
    }
    // m[0] = 0x01;
    // m[32] = 0x02;
    // m[64] = 0x04;
    // m[96] = 0x08;
    printf("--------------- masked message test ---------------------\n");
    unmask_bitstring(&m, 32);

    masked_poly r;
    poly_masked_frommsg(&r, &m);
    // printf("--------------- masked poly ---------------------\n");
    // print_masked_poly(&r);

    unsigned char rec_m[128];

    poly_masked_tomsg(&rec_m, &r);
    unmask_bitstring(&rec_m, 32);
    
    // -----------------------------------------------------------------------------
    // unsigned char unmasked_m[32] = {0x0f, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
    //                                 0x00, 0x04, 0x00, 0x04, 0x00, 0x1c, 0x00, 0x1c, 
    //                                 0x00, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
    //                                 0x00, 0x04, 0x00, 0x04, 0x00, 0x3c, 0x00, 0x3c};
    // printf("--------------- unmasked message test ---------------------\n");
    // print_bitstring(&unmasked_m);
    // poly unmasked_r;
    // poly_frommsg(&unmasked_r, &unmasked_m);
    // // printf("--------------- unmasked poly ---------------------\n");
    // // print_poly(&unmasked_r);
    // unsigned char unmasked_rec_m[32];
    // poly_tomsg(&unmasked_rec_m, &unmasked_r);
    // print_bitstring(&unmasked_rec_m);
}

void test_boolean_mask()
{
    unsigned char masked_m[33*4];
    unsigned char m[33] = { 0x0f, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
                            0x00, 0x04, 0x00, 0x04, 0x00, 0x1c, 0x00, 0x1c, 
                            0x00, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
                            0x00, 0x04, 0x00, 0x04, 0x00, 0x3c, 0x00, 0x3c, 0x11 };
    print_bitstring(&m, 33);
    random_boolean_mask(&masked_m, &m, 33);
    unmask_bitstring(&masked_m, 33);

}

void test_CBD()
{
    unsigned char seeds[32], masked_seeds[32 * 4];
    int nonce = 0, i;

    poly r;
    masked_poly masked_r;

    for (i = 0; i < 32; i++)
        seeds[i] = rand8();

    random_boolean_mask(masked_seeds, seeds,32);

    poly_sample(&r, seeds, nonce);

    print_poly(&r);

    poly_masked_sample(&masked_r, masked_seeds, nonce);

    print_masked_poly(&masked_r);
}

void test_enc()
{
    unsigned char pk[NEWHOPE_CPAPKE_PUBLICKEYBYTES], sk[NEWHOPE_CPAPKE_SECRETKEYBYTES];
    unsigned char c[NEWHOPE_CPAPKE_CIPHERTEXTBYTES], masked_c[NEWHOPE_CPAPKE_CIPHERTEXTBYTES];
    unsigned char m[32], rec_m[32], masked_m[32 * 4], masked_rec_m[32 * 4];
    unsigned char seed[32], masked_seed[32 * 4];
    masked_poly mskp;
    for (int i = 0; i < 32; i++)
    {
        m[i] = rand8();
        seed[i] = rand8();
    }
    printf("---------------- unmasked -------------\n");
    printf(" plaintext: ");
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
    printf(" recovered: ");
    print_bitstring(rec_m, 32);

    // masked_enc
    printf("---------------- masked -------------\n");
    printf(" plaintext: ");
    unmask_bitstring(masked_m, 32);
    cpapke_masked_enc(masked_c, masked_m, pk, masked_seed);
    printf(" ciphertext: ");
    print_bitstring(masked_c, 32);
    cpapke_masked_dec(masked_rec_m, masked_c, &mskp);
    printf(" recovered: ");
    unmask_bitstring(masked_rec_m, 32);



    

}

void test_frombytes()
{
    unsigned char a[NEWHOPE_CPAPKE_SECRETKEYBYTES], masked_a[NEWHOPE_CPAPKE_SECRETKEYBYTES * 4], rec_a[NEWHOPE_CPAPKE_SECRETKEYBYTES];
    poly r;
    masked_poly masked_r;
    for (int i = 0; i < NEWHOPE_CPAPKE_SECRETKEYBYTES; i++)
        a[i] = rand8();
    // printf("----------------------------------------\n");
    print_bitstring(a, NEWHOPE_CPAPKE_SECRETKEYBYTES);
    printf("----------------------------------------\n");
    random_boolean_mask(masked_a, a, NEWHOPE_CPAPKE_SECRETKEYBYTES);

    poly_frombytes(&r, a);
    poly_tobytes(rec_a, &r);
    print_bitstring(rec_a, NEWHOPE_CPAPKE_SECRETKEYBYTES);

    // print_poly(&r);
    // printf("----------------------------------------\n");
    poly_masked_frombytes(&masked_r, masked_a);
    // print_masked_poly(&masked_r);
}


void main()
{
    // test_msg();
    // test_boolean_mask();
    // test_CBD();
    test_enc();
    // test_frombytes();
}