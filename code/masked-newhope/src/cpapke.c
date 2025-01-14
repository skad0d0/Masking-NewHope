#include <stdio.h>
#include "api.h"
#include "poly.h"
#include "utils.h"
#include "randombytes.h"
#include "fips202.h"

/*************************************************
* Name:        encode_pk
* 
* Description: Serialize the public key as concatenation of the
*              serialization of the polynomial pk and the public seed
*              used to generete the polynomial a.
*
* Arguments:   unsigned char *r:          pointer to the output serialized public key
*              const poly *pk:            pointer to the input public-key polynomial
*              const unsigned char *seed: pointer to the input public seed
**************************************************/
static void encode_pk(unsigned char *r, const poly *pk, const unsigned char *seed)
{
  int i;
  poly_tobytes(r, pk);
  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    r[NEWHOPE_POLYBYTES+i] = seed[i];
}

/*************************************************
* Name:        decode_pk
* 
* Description: De-serialize the public key; inverse of encode_pk
*
* Arguments:   poly *pk:               pointer to output public-key polynomial
*              unsigned char *seed:    pointer to output public seed
*              const unsigned char *r: pointer to input byte array
**************************************************/
static void decode_pk(poly *pk, unsigned char *seed, const unsigned char *r)
{
  int i;
  poly_frombytes(pk, r);
  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    seed[i] = r[NEWHOPE_POLYBYTES+i];
}

/*************************************************
* Name:        encode_c
* 
* Description: Serialize the ciphertext as concatenation of the
*              serialization of the polynomial b and serialization
*              of the compressed polynomial v
*
* Arguments:   - unsigned char *r: pointer to the output serialized ciphertext
*              - const poly *b:    pointer to the input polynomial b
*              - const poly *v:    pointer to the input polynomial v
**************************************************/
static void encode_c(unsigned char *r, const poly *b, const poly *v)
{
  poly_tobytes(r,b);
  poly_compress(r+NEWHOPE_POLYBYTES,v);
}

/*************************************************
* Name:        decode_c
* 
* Description: de-serialize the ciphertext; inverse of encode_c
*
* Arguments:   - poly *b:                pointer to output polynomial b
*              - poly *v:                pointer to output polynomial v
*              - const unsigned char *r: pointer to input byte array
**************************************************/
void decode_c(poly *b, poly *v, const unsigned char *r)
{
  poly_frombytes(b, r);
  poly_decompress(v, r+NEWHOPE_POLYBYTES);
}

/*************************************************
* Name:        gen_a
* 
* Description: Deterministically generate public polynomial a from seed
*
* Arguments:   - poly *a:                   pointer to output polynomial a
*              - const unsigned char *seed: pointer to input seed
**************************************************/
static void gen_a(poly *a, const unsigned char *seed)
{
  poly_uniform(a,seed);
}


/*************************************************
* Name:        cpapke_keypair
* 
* Description: Generates public and private key 
*              for the CPA public-key encryption scheme underlying
*              the NewHope KEMs
*
* Arguments:   - unsigned char *pk: pointer to output public key
*              - unsigned char *sk: pointer to output private key
**************************************************/
void cpapke_keypair(unsigned char *pk,
                    unsigned char *sk)
{
  poly ahat, ehat, ahat_shat, bhat, shat;
  unsigned char z[2*NEWHOPE_SYMBYTES];
  unsigned char *publicseed = z;
  unsigned char *noiseseed = z+NEWHOPE_SYMBYTES;

  z[0] = 0x01;
  randombytes(z+1, NEWHOPE_SYMBYTES);
  shake256(z, 2*NEWHOPE_SYMBYTES, z, NEWHOPE_SYMBYTES + 1);

  gen_a(&ahat, publicseed);

  poly_sample(&shat, noiseseed, 0);
  poly_ntt(&shat);

  poly_sample(&ehat, noiseseed, 1);
  poly_ntt(&ehat);

  poly_mul_pointwise(&ahat_shat, &shat, &ahat);
  poly_add(&bhat, &ehat, &ahat_shat);

  poly_tobytes(sk, &shat);
  encode_pk(pk, &bhat, publicseed);
}

/*************************************************
* Name:        cpapke_enc
* 
* Description: Encryption function of
*              the CPA public-key encryption scheme underlying
*              the NewHope KEMs
*
* Arguments:   - unsigned char *c:          pointer to output ciphertext
*              - const unsigned char *m:    pointer to input message (of length NEWHOPE_SYMBYTES bytes)
*              - const unsigned char *pk:   pointer to input public key
*              - const unsigned char *coin: pointer to input random coins used as seed
*                                           to deterministically generate all randomness
**************************************************/
void cpapke_enc(unsigned char *c,
                const unsigned char *m,
                const unsigned char *pk,
                const unsigned char *coin)
{
  poly sprime, eprime, vprime, ahat, bhat, eprimeprime, uhat, v;
  unsigned char publicseed[NEWHOPE_SYMBYTES];

  poly_frommsg(&v, m); // need masking

  decode_pk(&bhat, publicseed, pk); // no masking
  gen_a(&ahat, publicseed);

  poly_sample(&sprime, coin, 0); // need masking
  poly_sample(&eprime, coin, 1);
  poly_sample(&eprimeprime, coin, 2);

  poly_ntt(&sprime); // need masking
  poly_ntt(&eprime);

  poly_mul_pointwise(&uhat, &ahat, &sprime);
  poly_add(&uhat, &uhat, &eprime);

  poly_mul_pointwise(&vprime, &bhat, &sprime);
  poly_invntt(&vprime);

  poly_add(&vprime, &vprime, &eprimeprime);
  poly_add(&vprime, &vprime, &v); // add message

  encode_c(c, &uhat, &vprime); // no masking 
}


/*************************************************
* Name:        cpapke_dec
* 
* Description: Decryption function of
*              the CPA public-key encryption scheme underlying
*              the NewHope KEMs
*
* Arguments:   - unsigned char *m:        pointer to output decrypted message
*              - const unsigned char *c:  pointer to input ciphertext
*              - const unsigned char *sk: pointer to input secret key
**************************************************/
void cpapke_dec(unsigned char *m,
                const unsigned char *c,
                const unsigned char *sk)
{
  poly vprime, uhat, tmp, shat;

  poly_frombytes(&shat, sk); // masking

  decode_c(&uhat, &vprime, c);
  poly_mul_pointwise(&tmp, &shat, &uhat);
  poly_invntt(&tmp);
  poly_sub(&tmp, &tmp, &vprime);
  poly_tomsg(m, &tmp); // masking
}

/**
 * @brief Encrypts a message using the masked CPA-secure public key encryption scheme.
 *
 * This function performs the encryption of a message using the masked CPA-secure public key encryption scheme.
 * It takes the masked message, public key, and masked coins as inputs, and produces the masked ciphertext.
 *
 * @param[out] masked_c       The output masked ciphertext.
 * @param[in]  masked_m       The input masked message.
 * @param[in]  pk             The public key.
 * @param[in]  masked_coins   The masked coins used for encryption.
 */
void cpapke_masked_enc(unsigned char *masked_c,
                 const unsigned char *masked_m,
                 const unsigned char *pk,
                 const unsigned char *masked_coins)
{
  poly bhat, ahat, uhat, vprime;
  masked_poly masked_v, masked_sprime, masked_eprime, masked_eprimeprime, masked_uhat, masked_vprime;
  unsigned char publicseed[NEWHOPE_SYMBYTES];

  poly_masked_frommsg(&masked_v, masked_m);
  decode_pk(&bhat, publicseed, pk);
  gen_a(&ahat, publicseed);

  poly_masked_sample(&masked_sprime, masked_coins, 0);
  poly_masked_sample(&masked_eprime, masked_coins, 1);
  poly_masked_sample(&masked_eprimeprime, masked_coins, 2);

  poly_masked_ntt(&masked_sprime);
  poly_masked_ntt(&masked_eprime);

  poly_halfmasked_mul_pointwise(&masked_uhat, &ahat, &masked_sprime);
  poly_masked_add(&masked_uhat, &masked_uhat, &masked_eprime);

  poly_halfmasked_mul_pointwise(&masked_vprime, &bhat, &masked_sprime);
  poly_masked_invntt(&masked_vprime);

  poly_masked_add(&masked_vprime, &masked_vprime, &masked_eprimeprime);
  poly_masked_add(&masked_vprime, &masked_vprime, &masked_v);

  unmask_poly(&masked_uhat, &uhat);
  unmask_poly(&masked_vprime, &vprime);

  encode_c(masked_c, &uhat, &vprime);
}

void cpapke_masked_enc_no_encode(Masked *mct,
                    const unsigned char *masked_m,
                    const unsigned char *pk,
                    const unsigned char *masked_coins)
{
  poly bhat, ahat, uhat, vprime;
  masked_poly masked_v, masked_sprime, masked_eprime, masked_eprimeprime, masked_uhat, masked_vprime;
  unsigned char publicseed[NEWHOPE_SYMBYTES];
  int i, k;
  

  poly_masked_frommsg(&masked_v, masked_m);
  decode_pk(&bhat, publicseed, pk);
  gen_a(&ahat, publicseed);

  poly_masked_sample(&masked_sprime, masked_coins, 0);
  poly_masked_sample(&masked_eprime, masked_coins, 1);
  poly_masked_sample(&masked_eprimeprime, masked_coins, 2);

  poly_masked_ntt(&masked_sprime);
  poly_masked_ntt(&masked_eprime);

  poly_halfmasked_mul_pointwise(&masked_uhat, &ahat, &masked_sprime);
  poly_masked_add(&masked_uhat, &masked_uhat, &masked_eprime);

  poly_halfmasked_mul_pointwise(&masked_vprime, &bhat, &masked_sprime);
  poly_masked_invntt(&masked_vprime);

  poly_masked_add(&masked_vprime, &masked_vprime, &masked_eprimeprime);
  poly_masked_add(&masked_vprime, &masked_vprime, &masked_v);

  for (i = 0; i < NEWHOPE_N; i++)
  {
    for (k = 0; k < NEWHOPE_MASKING_ORDER + 1; k++)
    {
      mct[    i].shares[k] = masked_uhat.poly_shares[k].coeffs[i];
      mct[512+i].shares[k] = masked_vprime.poly_shares[k].coeffs[i];
    }
  }
}



/**
 * @brief Decrypts a masked ciphertext using a masked secret key.
 *
 * This function performs the decryption of a masked ciphertext to produce
 * a masked message. It uses various polynomial operations to achieve this.
 *
 * @param[out] masked_m  Pointer to the output masked message.
 * @param[in]  masked_c  Pointer to the input masked ciphertext.
 * @param[in]  mskp      Pointer to the masked secret key polynomial.
 */
void cpapke_masked_dec(unsigned char *masked_m,
               const unsigned char *masked_c,
               masked_poly *mskp)
{
  masked_poly masked_tmp;
  poly uhat, vprime;
  
  // poly_masked_frombytes(&masked_shat, sk);

  decode_c(&uhat, &vprime, masked_c);
  
  poly_halfmasked_mul_pointwise(&masked_tmp, &uhat, mskp);
  poly_masked_invntt(&masked_tmp);
  poly_halfmasked_sub(&masked_tmp, &masked_tmp, &vprime);
  poly_masked_tomsg(masked_m, &masked_tmp);
}