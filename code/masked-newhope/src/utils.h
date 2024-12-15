#ifndef TEST_H
#define TEST_H
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "poly.h"
#include "gadgets.h"

void print_bitstring(unsigned char * bs);
void unmask_bitstring(unsigned char * bs);


void print_poly(poly* p);
// void compare_poly(poly* a, poly* b);
void print_masked_poly(masked_poly* mp);
void print_masked_poly_arith(masked_poly* x);
void print_masked_arith(Masked* x, int q);
void print_masked_bool(Masked* y);

#endif