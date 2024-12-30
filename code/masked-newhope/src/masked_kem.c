#include "masked_kem.h"

static unsigned mswitch(unsigned x, unsigned q_start, unsigned q_end){
  return (2*q_end*x+q_start)/(2*q_start);
}

static unsigned compress(unsigned x, unsigned q, unsigned d){
  return mswitch(x, q, 1<<d)%(1<<d);
}

int masked_crypto_kem_dec(unsigned char *ss, 
                    const unsigned char *ct, 
                            masked_poly *mskp,
                    const unsigned char *pk,
                    const unsigned char *pkh,
                    const unsigned char *masked_s)
{
    unsigned char buf[(2*NEWHOPE_SYMBYTES+1) * (NEWHOPE_MASKING_ORDER + 1)], m[32*(NEWHOPE_MASKING_ORDER + 1)];
    unsigned char masked_k_coin_d[3*NEWHOPE_SYMBYTES * (NEWHOPE_MASKING_ORDER + 1)];
    unsigned char masked_coin[32*(NEWHOPE_MASKING_ORDER + 1)], new_coin[32], masked_dprime[32*(NEWHOPE_MASKING_ORDER + 1)];
    unsigned char c[NEWHOPE_CPAPKE_CIPHERTEXTBYTES], d[32], masked_c[NEWHOPE_CPAPKE_CIPHERTEXTBYTES* (NEWHOPE_MASKING_ORDER + 1)];
    Masked mct[1024];
    poly uhat, vprime;
    uint16_t poly_ct[1024];
    int check;
    int i, k;

    for (i = 0; i < NEWHOPE_CPAPKE_CIPHERTEXTBYTES; i++) c[i] = ct[i];
    for (i = 0; i < 32; i++) d[i] = ct[i + NEWHOPE_CPAPKE_CIPHERTEXTBYTES];


    buf[0] = 0x08;
    for (k = 1; k < NEWHOPE_MASKING_ORDER + 1; k++)
        buf[0 + k*(2*NEWHOPE_SYMBYTES+1)] = 0;
    
    cpapke_masked_dec(m, c, mskp);

    for (i = 0; i < 32; i++)
    {
        buf[1 + i] = m[i];
        buf[1 + 32 + i] = pkh[i];
    }

    for (k = 1; k < NEWHOPE_MASKING_ORDER + 1; k++)
    {
        for (i = 0; i < 32; i++)
        {
            buf[1 + i + k * (2*NEWHOPE_SYMBYTES+1)] = m[i + 32*k];
            buf[1 + 32 + i + k * (2*NEWHOPE_SYMBYTES+1)] = 0;
        }
    }

    shake256_masked(masked_k_coin_d, 3*NEWHOPE_SYMBYTES, buf, 2*NEWHOPE_SYMBYTES+1);



    for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
    {
        for (i = 0; i < 32; i++)
        {
            masked_coin[i + 32*k] = masked_k_coin_d[i + 32 + k*96];
            masked_dprime[i + 32*k] = masked_k_coin_d[i + 64 + k*96];
        }
    }

    // cpapke_masked_enc_no_encode(m_uhat, m_vprime, m, pk, masked_coin);
    cpapke_masked_enc_no_encode(mct, m, pk, masked_coin);

    for (i = 0; i < 32; i++)
    {
        for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
            mct[1024 + i].shares[k] = masked_dprime[i + 32*k];
    }

    decode_c(&uhat, &vprime, c);

    for (i = 0; i < NEWHOPE_N; i++)
    {
        poly_ct[i      ] = compress(uhat.coeffs[i], NEWHOPE_Q, 8);
        poly_ct[i + 512] = compress(vprime.coeffs[i], NEWHOPE_Q, 3);
    }

    for (i = 0; i < 32; i++) 
        poly_ct[1024 + i] = d[i];
    

    // check_uhat = newhope_poly_comp_hybrid(m_uhat, uhat);
    // check_vprime = newhope_poly_comp_hybrid(m_vprime, vprime);
    // fail = 1 if either check_uhat or check_vprime is 0, fail = 0 if both are 1 
    // fail = !(check_uhat & check_vprime);
    check = newhope_poly_comp_hybrid(mct, poly_ct);
    /* overwrite coins in k_coin_d with H(c) */
    shake256(new_coin, 32, ct, NEWHOPE_CCAKEM_CIPHERTEXTBYTES);

    for (i = 0; i < 32; i++) masked_k_coin_d[32 + i] = new_coin[i];
    for (k = 1; k < NEWHOPE_MASKING_ORDER + 1; k++)
    {
        for (i = 0; i < 32; i++)
           masked_k_coin_d[32 + i + k*96] = 0;
    }

    /* Overwrite pre-k with s on re-encryption failure */
    if (!check)
    {
        for (i = 0; i < 32; i++)
        {
            for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
                masked_k_coin_d[i + k*96] = masked_s[i + 32*k];
        }
    }

    /* hash concatenation of pre-k and h(c) to k */
    shake256_masked(ss, 32, masked_k_coin_d, 96);
    return 0;
}