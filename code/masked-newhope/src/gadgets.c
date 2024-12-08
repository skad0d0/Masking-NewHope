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

