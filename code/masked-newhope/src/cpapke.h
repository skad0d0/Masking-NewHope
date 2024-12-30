#ifndef INDCPA_H
#define INDCPA_H

#include "gadgets.h"

void decode_c(poly *b, poly *v, const unsigned char *r);

void cpapke_keypair(unsigned char *pk, 
                    unsigned char *sk);

void cpapke_enc(unsigned char *c,
               const unsigned char *m,
               const unsigned char *pk,
               const unsigned char *coins);

void cpapke_dec(unsigned char *m,
               const unsigned char *c,
               const unsigned char *sk);

void cpapke_masked_enc_no_encode(Masked *m_uhat,
               Masked *m_vprime,
               const unsigned char *masked_m,
               const unsigned char *pk,
               const unsigned char *masked_coins);


void cpapke_masked_enc(unsigned char *masked_c,
               const unsigned char *masked_m,
               const unsigned char *pk,
               const unsigned char *masked_coins);

void cpapke_masked_dec(unsigned char *masked_m,
               const unsigned char *masked_c,
               masked_poly *mskp);

#endif
