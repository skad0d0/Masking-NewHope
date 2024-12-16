#include "poly.h"
#include "gadgets.h"
#include "randombytes.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

// test message encoding and decoding
void test_msg()
{
    unsigned char m[32*4];
    int i,k;
    for (i = 0; i < 32*4; i++)
    {
        m[i] = rand() & 0xff;
        // for (k = 0; k < 4; k++)
        //     m[i + 32*k] = i+k;
    }
    // m[0] = 0x01;
    // m[32] = 0x02;
    // m[64] = 0x04;
    // m[96] = 0x08;
    printf("--------------- masked message test ---------------------\n");
    unmask_bitstring(&m);

    masked_poly r;
    poly_masked_frommsg(&r, &m);
    // printf("--------------- masked poly ---------------------\n");
    // print_masked_poly(&r);

    unsigned char rec_m[128];

    poly_masked_tomsg(&rec_m, &r);
    unmask_bitstring(&rec_m);
    
    // -----------------------------------------------------------------------------
    // unsigned char unmasked_m[32] = {0x0f, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
    //                                 0x00, 0x04, 0x00, 0x04, 0x00, 0x1c, 0x00, 0x1c, 
    //                                 0x00, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
    //                                 0x00, 0x04, 0x00, 0x04, 0x00, 0x3c, 0x00, 0x3c};
    // printf("--------------- unmasked message test ---------------------\n");
    // print_bitstring(&unmasked_m);
    // poly unmasked_r;
    // poly_frommsg(&unmasked_r, &unmasked_m);
    // // printf("--------------- unmasked poly ---------------------\n");
    // // print_poly(&unmasked_r);
    // unsigned char unmasked_rec_m[32];
    // poly_tomsg(&unmasked_rec_m, &unmasked_r);
    // print_bitstring(&unmasked_rec_m);
}

void test_boolean_mask()
{
    unsigned char masked_m[128];
    unsigned char m[32] = { 0x0f, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
                            0x00, 0x04, 0x00, 0x04, 0x00, 0x1c, 0x00, 0x1c, 
                            0x00, 0x04, 0x00, 0x04, 0x00, 0x0c, 0x00, 0x0c, 
                            0x00, 0x04, 0x00, 0x04, 0x00, 0x3c, 0x00, 0x3c };
    print_bitstring(&m);
    random_boolean_mask(&masked_m, &m);
    unmask_bitstring(&masked_m);

}


void main()
{
    test_msg();
    // test_boolean_mask();
}