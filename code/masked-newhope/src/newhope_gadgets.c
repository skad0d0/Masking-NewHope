#include "gadgets.h"
#include <math.h>

// m is a boolean mask of the msg. 
void encode_message(const unsigned char *m, masked_poly* y)
{
    Masked t1,t2;
    for(int i=0; i < 32; ++i)
    {
        for(int j=0; j < 8; ++j)
        {
            for(int k=0; k < NEWHOPE_MASKING_ORDER+1; ++k) 
                t1.shares[k] = (m[i+k*32]>>j)&1; 
            
            convert_B2A(&t1, &t2);

            for(int k=0; k < NEWHOPE_MASKING_ORDER+1; ++k)
            {
                (y->poly_shares[k]).coeffs[i*8+j+0] = (t2.shares[k]*(NEWHOPE_Q/2)) % NEWHOPE_Q;
                (y->poly_shares[k]).coeffs[i*8+j+256] = (t2.shares[k]*(NEWHOPE_Q/2)) % NEWHOPE_Q;
            }
        }
    }
}

void modulus_switch(Masked* x, unsigned q, unsigned shift)
{
  /* 
   * Modulus switch between Z_q and Z_{2^shift} 
   * round((x<<shift)/q) = ((x<<(shift+1) + q1)//(2*q)
   * No overflow should appear for the values we use in the paper
   */
  int64_t temp;
  for(int i =0; i < NEWHOPE_MASKING_ORDER+1; ++i) {
    temp = (int64_t)(x->shares[i]) << (shift+1);
    temp = (temp+q)/(2*q);
    x->shares[i] = (int)temp&((1<<shift)-1);
  }
}


// unsigned switch_table[9] =  {6, 7, 7, 7, 8, 8, 8, 8, 8}; // Value of \ell in the paper

void newhope_decryption(Masked* x, Masked* b)
{
  // unsigned l = switch_table[NEWHOPE_MASKING_ORDER-1];
  modulus_switch(x, NEWHOPE_Q, 7);
  convert_2_l_to_1bit_bool(x, b);
}

void CBD(Masked* a, Masked* b, Masked* y)
{
  Masked h_a, h_b;
  Masked t1, t2;
  int eta = 8; /* y in [-8, 8] */
  
  for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j) t1.shares[j] = (a->shares[j])&1;
  convert_B2A(&t1, &h_a);
  for(int i=1; i < eta; ++i){
    for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j) t1.shares[j] = (a->shares[j] >> i)&1;
    convert_B2A(&t1, &t2);
    for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j) h_a.shares[j] = (h_a.shares[j] + t2.shares[j])%NEWHOPE_Q;
  }

  for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j) t1.shares[j] = (b->shares[j])&1;
  convert_B2A(&t1, &h_b);
  for(int i=1; i < eta; ++i){
    for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j) t1.shares[j] = (b->shares[j] >> i)&1;
    convert_B2A(&t1, &t2);
    for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j) h_b.shares[j] = (h_b.shares[j] + t2.shares[j])%NEWHOPE_Q;
  }

  for(int i =0; i < NEWHOPE_MASKING_ORDER+1; ++i){
    y->shares[i] = (h_a.shares[i] + NEWHOPE_Q - h_b.shares[i])%NEWHOPE_Q;
  }
}