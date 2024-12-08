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
                (y->poly_shares[k]).coeffs[i*8+j+0] = (t2.shares[k]*((NEWHOPE_Q+1)/2)) % NEWHOPE_Q;
                (y->poly_shares[k]).coeffs[i*8+j+512] = (t2.shares[k]*((NEWHOPE_Q+1)/2)) % NEWHOPE_Q;
            }
        }
    }
}