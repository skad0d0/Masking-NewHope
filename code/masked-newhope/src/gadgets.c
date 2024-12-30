#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "gadgets.h"

// linear refresh mask
void linear_arithmetic_refresh(Masked* x, unsigned q){
  uint16_t r[2];
  for(int i=0; i < (NEWHOPE_MASKING_ORDER-1); i+=2){
    rand_q(r);
    x->shares[i] = (x->shares[i] + r[0])%q;
    x->shares[NEWHOPE_MASKING_ORDER] = (x->shares[NEWHOPE_MASKING_ORDER] - r[0] + q)%q;

    x->shares[i+1] = (x->shares[i+1] + r[1])%q;
    x->shares[NEWHOPE_MASKING_ORDER] = (x->shares[NEWHOPE_MASKING_ORDER] - r[1] + q)%q;
  }
  #if NEWHOPE_MASKING_ORDER%2 == 1
  rand_q(r);
  x->shares[NEWHOPE_MASKING_ORDER-1] = (x->shares[NEWHOPE_MASKING_ORDER-1] + r[0])%q;
  x->shares[NEWHOPE_MASKING_ORDER] = (x->shares[NEWHOPE_MASKING_ORDER] - r[0] + q)%q;
  #endif
}

void linear_boolean_refresh(Masked* x, unsigned k){
  int r;
  for(int i=0; i< NEWHOPE_MASKING_ORDER; ++i){
    r = (int) rand32() & ((1<<k)-1);
    x->shares[i] = (x->shares[i] ^ r);
    x->shares[NEWHOPE_MASKING_ORDER] = (x->shares[NEWHOPE_MASKING_ORDER] ^ r);
  }
}

void boolean_refresh(Masked* x, unsigned k){
  int r;
  for(int i=0; i< NEWHOPE_MASKING_ORDER+1; ++i){
    for(int j=i+1; j < NEWHOPE_MASKING_ORDER+1; ++j){
      r = (int) rand32() & ((1<<k)-1);
      x->shares[i] = (x->shares[i] ^ r);
      x->shares[j] = (x->shares[j] ^ r);
    }
  }
}

void fill_masked_mod_q(Masked* x){
  uint16_t r[2];
  for(int i=0; i < NEWHOPE_MASKING_ORDER; i+=2){
    rand_q(r);
    x->shares[i] = r[0];
    x->shares[i+1] = r[1];
  }
  #if NEWHOPE_MASKING_ORDER%2 == 0
  rand_q(r);
  x->shares[NEWHOPE_MASKING_ORDER] = r[0];
  #endif
}

// convert boolean representation x to arithmetic representation y
void convert_B2A(Masked *x, Masked *y)
{
    Masked T[2];
    Masked T_p[2];
    int q = NEWHOPE_Q;

    for(int u=0; u < 2; ++u){
        T[u].shares[0] = u%q;
        for(int i=1; i < NEWHOPE_MASKING_ORDER+1; ++i) T[u].shares[i] = 0;
    }

    for(int i=0; i < NEWHOPE_MASKING_ORDER; ++i){
        for(int u=0; u < 2; ++u)
        {
            for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j)
                T_p[u].shares[j] = T[u^(x->shares[i])].shares[j]; 
        }
        for(int u=0; u < 2; ++u)
        {
            linear_arithmetic_refresh(&(T_p[u]), q); 
            for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j) 
                T[u].shares[j] = T_p[u].shares[j];
        }
    }

    for(int i=0; i < NEWHOPE_MASKING_ORDER+1; ++i) 
        y->shares[i] = T[x->shares[NEWHOPE_MASKING_ORDER]].shares[i]; 
    linear_arithmetic_refresh(y, q);
}

void convert_2_l_to_1bit_bool(Masked* x, Masked* b)
{
    uint64_t T1[NEWHOPE_MASKING_ORDER+1];
    uint64_t T2[NEWHOPE_MASKING_ORDER+1];
    uint64_t r;
    unsigned shift;

    T1[0] = 0xFFFFFFFF00000000LLU;
    T2[0] = 0x00000000FFFFFFFFLLU;
 

    for(int i=1; i < NEWHOPE_MASKING_ORDER+1; ++i) {
      T1[i] = 0;
      T2[i] = 0;
    }
    
    for(int i=0; i < NEWHOPE_MASKING_ORDER; ++i){
      for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j){
        shift = x->shares[i];
        if (shift%64 != 0){
          r = T1[j];
          T1[j] = (T1[j] >> (shift%64)) + (T2[j] << (64-(shift%64)));
          T2[j] = (T2[j] >> (shift%64)) + (r << (64-(shift%64)));
        }
        if (shift >= 64){
          r = T2[j];
          T2[j] = T1[j];
          T1[j] = r;
         }
      }
      for(int j=0; j < NEWHOPE_MASKING_ORDER; ++j){
        r = rand64();
        T1[j] ^= r;
        T1[NEWHOPE_MASKING_ORDER] ^= r;
        r = rand64();
        T2[j] ^= r;
        T2[NEWHOPE_MASKING_ORDER] ^= r;
      }
    }
    
    for(int i=0; i < NEWHOPE_MASKING_ORDER+1; ++i){
      shift = x->shares[NEWHOPE_MASKING_ORDER];
      if (shift < 64) b->shares[i] = (T1[i]>>(shift))&1;
      else            b->shares[i] = (T2[i]>>(shift-64))&1;
    }
    for(int j=0; j < NEWHOPE_MASKING_ORDER; ++j){
      r = rand32();
      b->shares[j] ^= r;
      b->shares[NEWHOPE_MASKING_ORDER] ^= r;
      b->shares[j] &= 1;
      b->shares[NEWHOPE_MASKING_ORDER] &= 1;
    }
  }

void random_boolean_mask(unsigned char *masked_m, unsigned char *m, int length)
{
  int i, k;
  unsigned char share;

  // generate random 8bits mask
  for (i = 0; i < length; i++)
  {
    share = m[i];
    for (k = 0; k < NEWHOPE_MASKING_ORDER; k++)
    {
      masked_m[i + length*k] = rand8();
      share ^= masked_m[i + length*k];
    }
    masked_m[i + length*NEWHOPE_MASKING_ORDER] = share;
  }

}
