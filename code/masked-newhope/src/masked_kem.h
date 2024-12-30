#ifndef MKEM_H
#define MKEM_H

#include "params.h"
#include "poly.h"
#include "cpapke.h"
#include "fips202.h"
// #include "gadgets.h"

int masked_crypto_kem_dec(unsigned char *ss, 
                    const unsigned char *ct, 
                            masked_poly *mskp,
                    const unsigned char *pk,
                    const unsigned char *pkh,
                    const unsigned char *masked_s);

#endif