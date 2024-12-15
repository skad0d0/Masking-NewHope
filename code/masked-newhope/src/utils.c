#include "utils.h"
#define DISPLAY 20
void print_masked_poly(masked_poly* mp){
  poly p;
  unmask_poly(mp,&p);
  print_poly(&p);
}

void print_poly(poly* p){
  for(int i=0; i < DISPLAY; ++i){
    printf("%i ", p->coeffs[i]); 
  }
  printf("\n");
}

void print_masked_arith(Masked* x, int q){
  printf(" (");
  int t=0;
  for(int i=0; i < NEWHOPE_MASKING_ORDER; ++i){
    printf("%i, ",x->shares[i]);
    t += x->shares[i];
  }
  t += x->shares[NEWHOPE_MASKING_ORDER];
  printf("%i) = %i = %i mod %u\n", x->shares[NEWHOPE_MASKING_ORDER], t, t%q, q);
}

void print_masked_poly_arith(masked_poly* x){
  int q = NEWHOPE_Q;
  int size = 10;
  Masked t;
  for(int i=0; i < size; ++i){
    for(int j=0; j < NEWHOPE_MASKING_ORDER+1; ++j){
      t.shares[j] = (x->poly_shares[j]).coeffs[i];
    }
    print_masked_arith(&t,q);
  }
}

void print_masked_bool(Masked* y){
  int t=0;
  printf(" (");
  for(int i=0; i < NEWHOPE_MASKING_ORDER; ++i){
    printf("%i, ",y->shares[i]);
    t ^= y->shares[i];
  }
  t ^= y->shares[NEWHOPE_MASKING_ORDER];
  printf("%i) = (", y->shares[NEWHOPE_MASKING_ORDER]);
  for(int i=0; i < NEWHOPE_MASKING_ORDER; ++i){
    printf("0x%X, ", y->shares[i]);
  }
  printf("0x%x) = %i\n", y->shares[NEWHOPE_MASKING_ORDER], t);
}

void print_bitstring(unsigned char * bs){
  for(int i=0; i < 32; ++i) printf("%X ",bs[i]);
  printf("\n");
}

void unmask_bitstring(unsigned char * bs){
  
  unsigned char unmasked_buf[32];

  for(int i=0; i < 32; ++i) unmasked_buf[i] = 0;
  for(int k=0; k < NEWHOPE_MASKING_ORDER+1; ++k){
    for(int i=0; i < 32; ++i) unmasked_buf[i] ^= bs[i+k*32];
  }
  print_bitstring(unmasked_buf);
}