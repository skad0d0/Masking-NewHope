#include "random.h"
#include "params.h"

#if RNG_MODE == 1
static unsigned x=123456789, y=362436069, z=521288629;
#endif

uint32_t rand32(){

#ifdef COUNT
  count_rand++;
#endif

#if RNG_MODE == 1
  unsigned t;

  x ^= x << 16;
  x ^= x >> 5;
  x ^= x << 1;

  t = x;
  x = y;
  y = z;
  z = t ^ x ^ y; 
  return z;
#elif RNG_MODE == 2
  return rand();
#elif RNG_MODE == 0
  return 0;
#endif
}

// v[0] = r mod q, v[1] = floor(r/q)
void rand_q(uint16_t v[2]){
  uint32_t r;
  do{
    r = rand32();
  } while(r > (28U * NEWHOPE_Q*NEWHOPE_Q)); /* rejection prob = 0.015 */
  r = r%(NEWHOPE_Q*NEWHOPE_Q);
  v[0] = r%(NEWHOPE_Q);
  v[1] = r/NEWHOPE_Q;
}