#ifndef FIPS202_H
#define FIPS202_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "params.h"


#define SHAKE128_RATE 168
#define SHAKE256_RATE 136

typedef struct {
  uint64_t s[25];
} keccak_state;

typedef struct {
	uint64_t s_masked[25 * (NEWHOPE_MASKING_ORDER + 1)];
} keccak_state_masked;


void shake128_absorb(uint64_t *s, const unsigned char *input, unsigned long long inputByteLen);
void shake128_squeezeblocks(unsigned char *output, unsigned long long nblocks, uint64_t *s);
void shake256(unsigned char *output, unsigned long long outputByteLen, const unsigned char *input, unsigned long long inputByteLen);

void shake256_absorb_masked(keccak_state_masked* state_masked, const unsigned char* in_masked, size_t inlen);
void shake256_squeezeblocks_masked(unsigned char* out_masked, size_t nblocks, keccak_state_masked* state_masked, size_t outlen);
void shake256_masked(unsigned char* out_masked, size_t outlen, const unsigned char* in_masked, size_t inlen);

#endif
