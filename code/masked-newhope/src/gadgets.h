#ifndef GADGETS_H
#define GADGETS_H

#include <stdint.h>
#include <stdio.h>

#include "params.h"
#include "poly.h"
#include "random.h"

typedef struct
{
     int shares[NEWHOPE_MASKING_ORDER+1];
} Masked;

void linear_arithmetic_refresh(Masked* x, unsigned q);
void linear_boolean_refresh(Masked* x, unsigned k);
void boolean_refresh(Masked* x, unsigned k);
void fill_masked_mod_q(Masked* x);
void convert_B2A(Masked *x, Masked *y);
unsigned char recombine_boolean_shares(Masked* x, int k);

void convert_2_l_to_1bit_bool(Masked* x, Masked* b);


void encode_message(const unsigned char *m, masked_poly* y);
void modulus_switch(Masked* x, unsigned q, unsigned shift);
void newhope_decryption(Masked* x, Masked* b);
void random_boolean_mask(unsigned char *masked_m, unsigned char *m, int length);
void CBD(Masked* a, Masked* b, Masked* y);
void sec_and(Masked* x, Masked* y, Masked* res, int k);

// high-order poly compare
int newhope_poly_comp_hybrid(Masked* mmasked_poly, uint16_t* ppoly);
int newhope_bool_comp(Masked *x, unsigned char *y);
#endif