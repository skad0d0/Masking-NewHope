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

void convert_B2A(Masked *x, Masked *y);

// unsigned char m[(32)*(NEWHOPE_MASKING_ORDER+1)]
void encode_message(const unsigned char *m, masked_poly* y);
#endif