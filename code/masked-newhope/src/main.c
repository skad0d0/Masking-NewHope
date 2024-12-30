#include "poly.h"
#include "gadgets.h"
#include "randombytes.h"
#include "cpapke.h"
#include "fips202.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <string.h>



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
    printf("----------------  masked  -------------\n");
    printf(" plaintext: ");
    unmask_bitstring(masked_m, 32);
    cpapke_masked_enc(masked_c, masked_m, pk, masked_seed);
    printf(" ciphertext: ");
    print_bitstring(masked_c, 32);
    cpapke_masked_dec(masked_rec_m, masked_c, &mskp);
    printf(" recovered: ");
    unmask_bitstring(masked_rec_m, 32);
}




void main()
{
    test_enc();

}