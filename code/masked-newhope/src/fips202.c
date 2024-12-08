/* Based on the public domain implementation in
 * crypto_hash/keccakc512/simple/ from http://bench.cr.yp.to/supercop.html
 * by Ronny Van Keer 
 * and the public domain "TweetFips202" implementation
 * from https://twitter.com/tweetfips202
 * by Gilles Van Assche, Daniel J. Bernstein, and Peter Schwabe */

#include <stdint.h>
#include <assert.h>
#include "fips202.h"

#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))
unsigned long rand32bits(void)
{
    unsigned long tmp_r;
    tmp_r = rand();
    tmp_r ^= rand() << 15;
    tmp_r ^= rand() << 30;
    return tmp_r;
}
//secMult from https://www.iacr.org/archive/ches2010/62250403/62250403.pdf
void secMult(uint64_t* c, uint64_t* a, uint64_t* b)
{
    unsigned int i, j, offset_i, offset_j;
    uint64_t r_ij[(NEWHOPE_MASKING_ORDER + 1) * (NEWHOPE_MASKING_ORDER + 1)];
    memset(r_ij, 0, (NEWHOPE_MASKING_ORDER + 1) * (NEWHOPE_MASKING_ORDER + 1) * 8);
    for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
    {
        offset_i = i * (NEWHOPE_MASKING_ORDER + 1);
        for (j = i + 1; j < (NEWHOPE_MASKING_ORDER + 1); j++)
        {
            offset_j = j * (NEWHOPE_MASKING_ORDER + 1);
            r_ij[j + offset_i] = ((uint64_t)(rand32bits()) << 32) + (uint64_t)(rand32bits());
            r_ij[i + offset_j] = (a[i] & b[j]);
            r_ij[i + offset_j] = r_ij[i + offset_j] ^ r_ij[j + offset_i];
            r_ij[i + offset_j] = r_ij[i + offset_j] ^ (a[j] & b[i]);
        }
    }
    for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
    {
        c[i] = a[i] & b[i];
        offset_i = i * (NEWHOPE_MASKING_ORDER + 1);
        for (j = 0; j < (NEWHOPE_MASKING_ORDER + 1); j++)
        {
            if (i != j)
                c[i] = c[i] ^ r_ij[j + offset_i];
        }
    }
    return;
}

void not_mult_xor(uint64_t* r, uint64_t* n, uint64_t* m, uint64_t* x)
{
    unsigned int i;
    //uint64_t tmp_share;
    r[0] = n[0] ^ 0xFFFFFFFFFFFFFFFF;
    for (i = 1; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        r[i] = n[i];
    secMult(r, r, m);
    for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
    {
        r[i] = r[i] ^ x[i];
    }
    return;
}

/*************************************************
* Name:        load64
* 
* Description: Load 8 bytes into uint64_t in little-endian order
*
* Arguments:   - const unsigned char *x: pointer to input byte array
*
* Returns the loaded 64-bit unsigned integer
**************************************************/
static uint64_t load64(const unsigned char *x)
{
  unsigned long long r = 0, i;

  for (i = 0; i < 8; ++i) {
    r |= (unsigned long long)x[i] << 8 * i;
  }
  return r;
}

/*************************************************
* Name:        store64
* 
* Description: Store a 64-bit integer to a byte array in little-endian order
*
* Arguments:   - uint8_t *x: pointer to the output byte array
*              - uint64_t u: input 64-bit unsigned integer
**************************************************/
static void store64(uint8_t *x, uint64_t u)
{
  unsigned int i;

  for(i=0; i<8; ++i) {
    x[i] = u;
    u >>= 8;
  }
}

/* Keccak round constants */
static const uint64_t KeccakF_RoundConstants[NROUNDS] = 
{
    (uint64_t)0x0000000000000001ULL,
    (uint64_t)0x0000000000008082ULL,
    (uint64_t)0x800000000000808aULL,
    (uint64_t)0x8000000080008000ULL,
    (uint64_t)0x000000000000808bULL,
    (uint64_t)0x0000000080000001ULL,
    (uint64_t)0x8000000080008081ULL,
    (uint64_t)0x8000000000008009ULL,
    (uint64_t)0x000000000000008aULL,
    (uint64_t)0x0000000000000088ULL,
    (uint64_t)0x0000000080008009ULL,
    (uint64_t)0x000000008000000aULL,
    (uint64_t)0x000000008000808bULL,
    (uint64_t)0x800000000000008bULL,
    (uint64_t)0x8000000000008089ULL,
    (uint64_t)0x8000000000008003ULL,
    (uint64_t)0x8000000000008002ULL,
    (uint64_t)0x8000000000000080ULL,
    (uint64_t)0x000000000000800aULL,
    (uint64_t)0x800000008000000aULL,
    (uint64_t)0x8000000080008081ULL,
    (uint64_t)0x8000000000008080ULL,
    (uint64_t)0x0000000080000001ULL,
    (uint64_t)0x8000000080008008ULL
};

void KeccakF1600_StatePermute_masked(uint64_t state_masked[25 * (NEWHOPE_MASKING_ORDER + 1)])
{
    int round, i;

    uint64_t Aba[(NEWHOPE_MASKING_ORDER + 1)], Abe[(NEWHOPE_MASKING_ORDER + 1)], Abi[(NEWHOPE_MASKING_ORDER + 1)], Abo[(NEWHOPE_MASKING_ORDER + 1)], Abu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Aga[(NEWHOPE_MASKING_ORDER + 1)], Age[(NEWHOPE_MASKING_ORDER + 1)], Agi[(NEWHOPE_MASKING_ORDER + 1)], Ago[(NEWHOPE_MASKING_ORDER + 1)], Agu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Aka[(NEWHOPE_MASKING_ORDER + 1)], Ake[(NEWHOPE_MASKING_ORDER + 1)], Aki[(NEWHOPE_MASKING_ORDER + 1)], Ako[(NEWHOPE_MASKING_ORDER + 1)], Aku[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Ama[(NEWHOPE_MASKING_ORDER + 1)], Ame[(NEWHOPE_MASKING_ORDER + 1)], Ami[(NEWHOPE_MASKING_ORDER + 1)], Amo[(NEWHOPE_MASKING_ORDER + 1)], Amu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Asa[(NEWHOPE_MASKING_ORDER + 1)], Ase[(NEWHOPE_MASKING_ORDER + 1)], Asi[(NEWHOPE_MASKING_ORDER + 1)], Aso[(NEWHOPE_MASKING_ORDER + 1)], Asu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t BCa[(NEWHOPE_MASKING_ORDER + 1)], BCe[(NEWHOPE_MASKING_ORDER + 1)], BCi[(NEWHOPE_MASKING_ORDER + 1)], BCo[(NEWHOPE_MASKING_ORDER + 1)], BCu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Da[(NEWHOPE_MASKING_ORDER + 1)],  De[(NEWHOPE_MASKING_ORDER + 1)],  Di[(NEWHOPE_MASKING_ORDER + 1)],  Do[(NEWHOPE_MASKING_ORDER + 1)],  Du[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Eba[(NEWHOPE_MASKING_ORDER + 1)], Ebe[(NEWHOPE_MASKING_ORDER + 1)], Ebi[(NEWHOPE_MASKING_ORDER + 1)], Ebo[(NEWHOPE_MASKING_ORDER + 1)], Ebu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Ega[(NEWHOPE_MASKING_ORDER + 1)], Ege[(NEWHOPE_MASKING_ORDER + 1)], Egi[(NEWHOPE_MASKING_ORDER + 1)], Ego[(NEWHOPE_MASKING_ORDER + 1)], Egu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Eka[(NEWHOPE_MASKING_ORDER + 1)], Eke[(NEWHOPE_MASKING_ORDER + 1)], Eki[(NEWHOPE_MASKING_ORDER + 1)], Eko[(NEWHOPE_MASKING_ORDER + 1)], Eku[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Ema[(NEWHOPE_MASKING_ORDER + 1)], Eme[(NEWHOPE_MASKING_ORDER + 1)], Emi[(NEWHOPE_MASKING_ORDER + 1)], Emo[(NEWHOPE_MASKING_ORDER + 1)], Emu[(NEWHOPE_MASKING_ORDER + 1)];
    uint64_t Esa[(NEWHOPE_MASKING_ORDER + 1)], Ese[(NEWHOPE_MASKING_ORDER + 1)], Esi[(NEWHOPE_MASKING_ORDER + 1)], Eso[(NEWHOPE_MASKING_ORDER + 1)], Esu[(NEWHOPE_MASKING_ORDER + 1)];

    //copyFromState(A, state)
    for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
    {
        Aba[i] = state_masked[ 0+(i*25)];
        Abe[i] = state_masked[ 1+(i*25)];
        Abi[i] = state_masked[ 2+(i*25)];
        Abo[i] = state_masked[ 3+(i*25)];
        Abu[i] = state_masked[ 4+(i*25)];
        Aga[i] = state_masked[ 5+(i*25)];
        Age[i] = state_masked[ 6+(i*25)];
        Agi[i] = state_masked[ 7+(i*25)];
        Ago[i] = state_masked[ 8+(i*25)];
        Agu[i] = state_masked[ 9+(i*25)];
        Aka[i] = state_masked[10+(i*25)];
        Ake[i] = state_masked[11+(i*25)];
        Aki[i] = state_masked[12+(i*25)];
        Ako[i] = state_masked[13+(i*25)];
        Aku[i] = state_masked[14+(i*25)];
        Ama[i] = state_masked[15+(i*25)];
        Ame[i] = state_masked[16+(i*25)];
        Ami[i] = state_masked[17+(i*25)];
        Amo[i] = state_masked[18+(i*25)];
        Amu[i] = state_masked[19+(i*25)];
        Asa[i] = state_masked[20+(i*25)];
        Ase[i] = state_masked[21+(i*25)];
        Asi[i] = state_masked[22+(i*25)];
        Aso[i] = state_masked[23+(i*25)];
        Asu[i] = state_masked[24+(i*25)];
    }
    for (round = 0; round < NROUNDS; round += 2)
    {
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            //    prepareTheta
            BCa[i] = Aba[i] ^ Aga[i] ^ Aka[i] ^ Ama[i] ^ Asa[i];
            BCe[i] = Abe[i] ^ Age[i] ^ Ake[i] ^ Ame[i] ^ Ase[i];
            BCi[i] = Abi[i] ^ Agi[i] ^ Aki[i] ^ Ami[i] ^ Asi[i];
            BCo[i] = Abo[i] ^ Ago[i] ^ Ako[i] ^ Amo[i] ^ Aso[i];
            BCu[i] = Abu[i] ^ Agu[i] ^ Aku[i] ^ Amu[i] ^ Asu[i];

            //thetaRhoPiChiIotaPrepareTheta(round  , A, E)
            Da[i] = BCu[i] ^ ROL(BCe[i], 1);
            De[i] = BCa[i] ^ ROL(BCi[i], 1);
            Di[i] = BCe[i] ^ ROL(BCo[i], 1);
            Do[i] = BCi[i] ^ ROL(BCu[i], 1);
            Du[i] = BCo[i] ^ ROL(BCa[i], 1);

            Aba[i] ^= Da[i];
            BCa[i] = Aba[i];
            Age[i] ^= De[i];
            BCe[i] = ROL(Age[i], 44);
            Aki[i] ^= Di[i];
            BCi[i] = ROL(Aki[i], 43);
            Amo[i] ^= Do[i];
            BCo[i] = ROL(Amo[i], 21);
            Asu[i] ^= Du[i];
            BCu[i] = ROL(Asu[i], 14);
        }
        not_mult_xor(Eba, BCe, BCi, BCa);
        Eba[0] ^= (uint64_t)KeccakF_RoundConstants[round];
        not_mult_xor(Ebe, BCi, BCo, BCe);
        not_mult_xor(Ebi, BCo, BCu, BCi);
        not_mult_xor(Ebo, BCu, BCa, BCo);
        not_mult_xor(Ebu, BCa, BCe, BCu);
        //Eba[i] = BCa ^ ((~BCe) & BCi);
        //Eba[i] ^= (uint64_t)KeccakF_RoundConstants[round];
        //Ebe[i] = BCe ^ ((~BCi) & BCo);
        //Ebi[i] = BCi ^ ((~BCo) & BCu);
        //Ebo[i] = BCo ^ ((~BCu) & BCa);
        //Ebu[i] = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Abo[i] ^= Do[i];
            BCa[i] = ROL(Abo[i], 28);
            Agu[i] ^= Du[i];
            BCe[i] = ROL(Agu[i], 20);
            Aka[i] ^= Da[i];
            BCi[i] = ROL(Aka[i], 3);
            Ame[i] ^= De[i];
            BCo[i] = ROL(Ame[i], 45);
            Asi[i] ^= Di[i];
            BCu[i] = ROL(Asi[i], 61);
        }
        not_mult_xor(Ega, BCe, BCi, BCa);
        not_mult_xor(Ege, BCi, BCo, BCe);
        not_mult_xor(Egi, BCo, BCu, BCi);
        not_mult_xor(Ego, BCu, BCa, BCo);
        not_mult_xor(Egu, BCa, BCe, BCu);
        //Ega = BCa ^ ((~BCe) & BCi);
        //Ege = BCe ^ ((~BCi) & BCo);
        //Egi = BCi ^ ((~BCo) & BCu);
        //Ego = BCo ^ ((~BCu) & BCa);
        //Egu = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Abe[i] ^= De[i];
            BCa[i] = ROL(Abe[i], 1);
            Agi[i] ^= Di[i];
            BCe[i] = ROL(Agi[i], 6);
            Ako[i] ^= Do[i];
            BCi[i] = ROL(Ako[i], 25);
            Amu[i] ^= Du[i];
            BCo[i] = ROL(Amu[i], 8);
            Asa[i] ^= Da[i];
            BCu[i] = ROL(Asa[i], 18);
        }
        not_mult_xor(Eka, BCe, BCi, BCa);
        not_mult_xor(Eke, BCi, BCo, BCe);
        not_mult_xor(Eki, BCo, BCu, BCi);
        not_mult_xor(Eko, BCu, BCa, BCo);
        not_mult_xor(Eku, BCa, BCe, BCu);
        //Eka = BCa ^ ((~BCe) & BCi);
        //Eke = BCe ^ ((~BCi) & BCo);
        //Eki = BCi ^ ((~BCo) & BCu);
        //Eko = BCo ^ ((~BCu) & BCa);
        //Eku = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Abu[i] ^= Du[i];
            BCa[i] = ROL(Abu[i], 27);
            Aga[i] ^= Da[i];
            BCe[i] = ROL(Aga[i], 36);
            Ake[i] ^= De[i];
            BCi[i] = ROL(Ake[i], 10);
            Ami[i] ^= Di[i];
            BCo[i] = ROL(Ami[i], 15);
            Aso[i] ^= Do[i];
            BCu[i] = ROL(Aso[i], 56);
        }
        not_mult_xor(Ema, BCe, BCi, BCa);
        not_mult_xor(Eme, BCi, BCo, BCe);
        not_mult_xor(Emi, BCo, BCu, BCi);
        not_mult_xor(Emo, BCu, BCa, BCo);
        not_mult_xor(Emu, BCa, BCe, BCu);
        //Ema = BCa ^ ((~BCe) & BCi);
        //Eme = BCe ^ ((~BCi) & BCo);
        //Emi = BCi ^ ((~BCo) & BCu);
        //Emo = BCo ^ ((~BCu) & BCa);
        //Emu = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Abi[i] ^= Di[i];
            BCa[i] = ROL(Abi[i], 62);
            Ago[i] ^= Do[i];
            BCe[i] = ROL(Ago[i], 55);
            Aku[i] ^= Du[i];
            BCi[i] = ROL(Aku[i], 39);
            Ama[i] ^= Da[i];
            BCo[i] = ROL(Ama[i], 41);
            Ase[i] ^= De[i];
            BCu[i] = ROL(Ase[i], 2);
        }
        not_mult_xor(Esa, BCe, BCi, BCa);
        not_mult_xor(Ese, BCi, BCo, BCe);
        not_mult_xor(Esi, BCo, BCu, BCi);
        not_mult_xor(Eso, BCu, BCa, BCo);
        not_mult_xor(Esu, BCa, BCe, BCu);
        //Esa = BCa ^ ((~BCe) & BCi);
        //Ese = BCe ^ ((~BCi) & BCo);
        //Esi = BCi ^ ((~BCo) & BCu);
        //Eso = BCo ^ ((~BCu) & BCa);
        //Esu = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            //    prepareTheta
            BCa[i] = Eba[i] ^ Ega[i] ^ Eka[i] ^ Ema[i] ^ Esa[i];
            BCe[i] = Ebe[i] ^ Ege[i] ^ Eke[i] ^ Eme[i] ^ Ese[i];
            BCi[i] = Ebi[i] ^ Egi[i] ^ Eki[i] ^ Emi[i] ^ Esi[i];
            BCo[i] = Ebo[i] ^ Ego[i] ^ Eko[i] ^ Emo[i] ^ Eso[i];
            BCu[i] = Ebu[i] ^ Egu[i] ^ Eku[i] ^ Emu[i] ^ Esu[i];

            //thetaRhoPiChiIotaPrepareTheta(round+1, E, A)
            Da[i] = BCu[i] ^ ROL(BCe[i], 1);
            De[i] = BCa[i] ^ ROL(BCi[i], 1);
            Di[i] = BCe[i] ^ ROL(BCo[i], 1);
            Do[i] = BCi[i] ^ ROL(BCu[i], 1);
            Du[i] = BCo[i] ^ ROL(BCa[i], 1);

            Eba[i] ^= Da[i];
            BCa[i] = Eba[i];
            Ege[i] ^= De[i];
            BCe[i] = ROL(Ege[i], 44);
            Eki[i] ^= Di[i];
            BCi[i] = ROL(Eki[i], 43);
            Emo[i] ^= Do[i];
            BCo[i] = ROL(Emo[i], 21);
            Esu[i] ^= Du[i];
            BCu[i] = ROL(Esu[i], 14);
        }
        not_mult_xor(Aba, BCe, BCi, BCa);
        Aba[0] ^= (uint64_t)KeccakF_RoundConstants[round + 1];
        not_mult_xor(Abe, BCi, BCo, BCe);
        not_mult_xor(Abi, BCo, BCu, BCi);
        not_mult_xor(Abo, BCu, BCa, BCo);
        not_mult_xor(Abu, BCa, BCe, BCu);
        //Aba = BCa ^ ((~BCe) & BCi);
        //Aba ^= (uint64_t)KeccakF_RoundConstants[round + 1];
        //Abe = BCe ^ ((~BCi) & BCo);
        //Abi = BCi ^ ((~BCo) & BCu);
        //Abo = BCo ^ ((~BCu) & BCa);
        //Abu = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Ebo[i] ^= Do[i];
            BCa[i] = ROL(Ebo[i], 28);
            Egu[i] ^= Du[i];
            BCe[i] = ROL(Egu[i], 20);
            Eka[i] ^= Da[i];
            BCi[i] = ROL(Eka[i], 3);
            Eme[i] ^= De[i];
            BCo[i] = ROL(Eme[i], 45);
            Esi[i] ^= Di[i];
            BCu[i] = ROL(Esi[i], 61);
        }
        not_mult_xor(Aga, BCe, BCi, BCa);
        not_mult_xor(Age, BCi, BCo, BCe);
        not_mult_xor(Agi, BCo, BCu, BCi);
        not_mult_xor(Ago, BCu, BCa, BCo);
        not_mult_xor(Agu, BCa, BCe, BCu);
        //Aga = BCa ^ ((~BCe) & BCi);
        //Age = BCe ^ ((~BCi) & BCo);
        //Agi = BCi ^ ((~BCo) & BCu);
        //Ago = BCo ^ ((~BCu) & BCa);
        //Agu = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Ebe[i] ^= De[i];
            BCa[i] = ROL(Ebe[i], 1);
            Egi[i] ^= Di[i];
            BCe[i] = ROL(Egi[i], 6);
            Eko[i] ^= Do[i];
            BCi[i] = ROL(Eko[i], 25);
            Emu[i] ^= Du[i];
            BCo[i] = ROL(Emu[i], 8);
            Esa[i] ^= Da[i];
            BCu[i] = ROL(Esa[i], 18);
        }
        not_mult_xor(Aka, BCe, BCi, BCa);
        not_mult_xor(Ake, BCi, BCo, BCe);
        not_mult_xor(Aki, BCo, BCu, BCi);
        not_mult_xor(Ako, BCu, BCa, BCo);
        not_mult_xor(Aku, BCa, BCe, BCu);
        //Aka = BCa ^ ((~BCe) & BCi);
        //Ake = BCe ^ ((~BCi) & BCo);
        //Aki = BCi ^ ((~BCo) & BCu);
        //Ako = BCo ^ ((~BCu) & BCa);
        //Aku = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Ebu[i] ^= Du[i];
            BCa[i] = ROL(Ebu[i], 27);
            Ega[i] ^= Da[i];
            BCe[i] = ROL(Ega[i], 36);
            Eke[i] ^= De[i];
            BCi[i] = ROL(Eke[i], 10);
            Emi[i] ^= Di[i];
            BCo[i] = ROL(Emi[i], 15);
            Eso[i] ^= Do[i];
            BCu[i] = ROL(Eso[i], 56);
        }
        not_mult_xor(Ama, BCe, BCi, BCa);
        not_mult_xor(Ame, BCi, BCo, BCe);
        not_mult_xor(Ami, BCo, BCu, BCi);
        not_mult_xor(Amo, BCu, BCa, BCo);
        not_mult_xor(Amu, BCa, BCe, BCu);
        //Ama = BCa ^ ((~BCe) & BCi);
        //Ame = BCe ^ ((~BCi) & BCo);
        //Ami = BCi ^ ((~BCo) & BCu);
        //Amo = BCo ^ ((~BCu) & BCa);
        //Amu = BCu ^ ((~BCa) & BCe);
        for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
        {
            Ebi[i] ^= Di[i];
            BCa[i] = ROL(Ebi[i], 62);
            Ego[i] ^= Do[i];
            BCe[i] = ROL(Ego[i], 55);
            Eku[i] ^= Du[i];
            BCi[i] = ROL(Eku[i], 39);
            Ema[i] ^= Da[i];
            BCo[i] = ROL(Ema[i], 41);
            Ese[i] ^= De[i];
            BCu[i] = ROL(Ese[i], 2);
        }
        not_mult_xor(Asa, BCe, BCi, BCa);
        not_mult_xor(Ase, BCi, BCo, BCe);
        not_mult_xor(Asi, BCo, BCu, BCi);
        not_mult_xor(Aso, BCu, BCa, BCo);
        not_mult_xor(Asu, BCa, BCe, BCu);
        //Asa = BCa ^ ((~BCe) & BCi);
        //Ase = BCe ^ ((~BCi) & BCo);
        //Asi = BCi ^ ((~BCo) & BCu);
        //Aso = BCo ^ ((~BCu) & BCa);
        //Asu = BCu ^ ((~BCa) & BCe);
    }

    for (i = 0; i < (NEWHOPE_MASKING_ORDER + 1); i++)
    {
        //copyToState(state, A)
        state_masked[ 0+(i*25)] = Aba[i];
        state_masked[ 1+(i*25)] = Abe[i];
        state_masked[ 2+(i*25)] = Abi[i];
        state_masked[ 3+(i*25)] = Abo[i];
        state_masked[ 4+(i*25)] = Abu[i];
        state_masked[ 5+(i*25)] = Aga[i];
        state_masked[ 6+(i*25)] = Age[i];
        state_masked[ 7+(i*25)] = Agi[i];
        state_masked[ 8+(i*25)] = Ago[i];
        state_masked[ 9+(i*25)] = Agu[i];
        state_masked[10+(i*25)] = Aka[i];
        state_masked[11+(i*25)] = Ake[i];
        state_masked[12+(i*25)] = Aki[i];
        state_masked[13+(i*25)] = Ako[i];
        state_masked[14+(i*25)] = Aku[i];
        state_masked[15+(i*25)] = Ama[i];
        state_masked[16+(i*25)] = Ame[i];
        state_masked[17+(i*25)] = Ami[i];
        state_masked[18+(i*25)] = Amo[i];
        state_masked[19+(i*25)] = Amu[i];
        state_masked[20+(i*25)] = Asa[i];
        state_masked[21+(i*25)] = Ase[i];
        state_masked[22+(i*25)] = Asi[i];
        state_masked[23+(i*25)] = Aso[i];
        state_masked[24+(i*25)] = Asu[i];
    }
	return;
}



/*************************************************
* Name:        KeccakF1600_StatePermute
* 
* Description: The Keccak F1600 Permutation
*
* Arguments:   - uint64_t * state: pointer to in/output Keccak state
**************************************************/
void KeccakF1600_StatePermute(uint64_t * state)
{
  int round;

        uint64_t Aba, Abe, Abi, Abo, Abu;
        uint64_t Aga, Age, Agi, Ago, Agu;
        uint64_t Aka, Ake, Aki, Ako, Aku;
        uint64_t Ama, Ame, Ami, Amo, Amu;
        uint64_t Asa, Ase, Asi, Aso, Asu;
        uint64_t BCa, BCe, BCi, BCo, BCu;
        uint64_t Da, De, Di, Do, Du;
        uint64_t Eba, Ebe, Ebi, Ebo, Ebu;
        uint64_t Ega, Ege, Egi, Ego, Egu;
        uint64_t Eka, Eke, Eki, Eko, Eku;
        uint64_t Ema, Eme, Emi, Emo, Emu;
        uint64_t Esa, Ese, Esi, Eso, Esu;

        //copyFromState(A, state)
        Aba = state[ 0];
        Abe = state[ 1];
        Abi = state[ 2];
        Abo = state[ 3];
        Abu = state[ 4];
        Aga = state[ 5];
        Age = state[ 6];
        Agi = state[ 7];
        Ago = state[ 8];
        Agu = state[ 9];
        Aka = state[10];
        Ake = state[11];
        Aki = state[12];
        Ako = state[13];
        Aku = state[14];
        Ama = state[15];
        Ame = state[16];
        Ami = state[17];
        Amo = state[18];
        Amu = state[19];
        Asa = state[20];
        Ase = state[21];
        Asi = state[22];
        Aso = state[23];
        Asu = state[24];

        for( round = 0; round < NROUNDS; round += 2 )
        {
            //    prepareTheta
            BCa = Aba^Aga^Aka^Ama^Asa;
            BCe = Abe^Age^Ake^Ame^Ase;
            BCi = Abi^Agi^Aki^Ami^Asi;
            BCo = Abo^Ago^Ako^Amo^Aso;
            BCu = Abu^Agu^Aku^Amu^Asu;

            //thetaRhoPiChiIotaPrepareTheta(round  , A, E)
            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Aba ^= Da;
            BCa = Aba;
            Age ^= De;
            BCe = ROL(Age, 44);
            Aki ^= Di;
            BCi = ROL(Aki, 43);
            Amo ^= Do;
            BCo = ROL(Amo, 21);
            Asu ^= Du;
            BCu = ROL(Asu, 14);
            Eba =   BCa ^((~BCe)&  BCi );
            Eba ^= (uint64_t)KeccakF_RoundConstants[round];
            Ebe =   BCe ^((~BCi)&  BCo );
            Ebi =   BCi ^((~BCo)&  BCu );
            Ebo =   BCo ^((~BCu)&  BCa );
            Ebu =   BCu ^((~BCa)&  BCe );

            Abo ^= Do;
            BCa = ROL(Abo, 28);
            Agu ^= Du;
            BCe = ROL(Agu, 20);
            Aka ^= Da;
            BCi = ROL(Aka,  3);
            Ame ^= De;
            BCo = ROL(Ame, 45);
            Asi ^= Di;
            BCu = ROL(Asi, 61);
            Ega =   BCa ^((~BCe)&  BCi );
            Ege =   BCe ^((~BCi)&  BCo );
            Egi =   BCi ^((~BCo)&  BCu );
            Ego =   BCo ^((~BCu)&  BCa );
            Egu =   BCu ^((~BCa)&  BCe );

            Abe ^= De;
            BCa = ROL(Abe,  1);
            Agi ^= Di;
            BCe = ROL(Agi,  6);
            Ako ^= Do;
            BCi = ROL(Ako, 25);
            Amu ^= Du;
            BCo = ROL(Amu,  8);
            Asa ^= Da;
            BCu = ROL(Asa, 18);
            Eka =   BCa ^((~BCe)&  BCi );
            Eke =   BCe ^((~BCi)&  BCo );
            Eki =   BCi ^((~BCo)&  BCu );
            Eko =   BCo ^((~BCu)&  BCa );
            Eku =   BCu ^((~BCa)&  BCe );

            Abu ^= Du;
            BCa = ROL(Abu, 27);
            Aga ^= Da;
            BCe = ROL(Aga, 36);
            Ake ^= De;
            BCi = ROL(Ake, 10);
            Ami ^= Di;
            BCo = ROL(Ami, 15);
            Aso ^= Do;
            BCu = ROL(Aso, 56);
            Ema =   BCa ^((~BCe)&  BCi );
            Eme =   BCe ^((~BCi)&  BCo );
            Emi =   BCi ^((~BCo)&  BCu );
            Emo =   BCo ^((~BCu)&  BCa );
            Emu =   BCu ^((~BCa)&  BCe );

            Abi ^= Di;
            BCa = ROL(Abi, 62);
            Ago ^= Do;
            BCe = ROL(Ago, 55);
            Aku ^= Du;
            BCi = ROL(Aku, 39);
            Ama ^= Da;
            BCo = ROL(Ama, 41);
            Ase ^= De;
            BCu = ROL(Ase,  2);
            Esa =   BCa ^((~BCe)&  BCi );
            Ese =   BCe ^((~BCi)&  BCo );
            Esi =   BCi ^((~BCo)&  BCu );
            Eso =   BCo ^((~BCu)&  BCa );
            Esu =   BCu ^((~BCa)&  BCe );

            //    prepareTheta
            BCa = Eba^Ega^Eka^Ema^Esa;
            BCe = Ebe^Ege^Eke^Eme^Ese;
            BCi = Ebi^Egi^Eki^Emi^Esi;
            BCo = Ebo^Ego^Eko^Emo^Eso;
            BCu = Ebu^Egu^Eku^Emu^Esu;

            //thetaRhoPiChiIotaPrepareTheta(round+1, E, A)
            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Eba ^= Da;
            BCa = Eba;
            Ege ^= De;
            BCe = ROL(Ege, 44);
            Eki ^= Di;
            BCi = ROL(Eki, 43);
            Emo ^= Do;
            BCo = ROL(Emo, 21);
            Esu ^= Du;
            BCu = ROL(Esu, 14);
            Aba =   BCa ^((~BCe)&  BCi );
            Aba ^= (uint64_t)KeccakF_RoundConstants[round+1];
            Abe =   BCe ^((~BCi)&  BCo );
            Abi =   BCi ^((~BCo)&  BCu );
            Abo =   BCo ^((~BCu)&  BCa );
            Abu =   BCu ^((~BCa)&  BCe );

            Ebo ^= Do;
            BCa = ROL(Ebo, 28);
            Egu ^= Du;
            BCe = ROL(Egu, 20);
            Eka ^= Da;
            BCi = ROL(Eka, 3);
            Eme ^= De;
            BCo = ROL(Eme, 45);
            Esi ^= Di;
            BCu = ROL(Esi, 61);
            Aga =   BCa ^((~BCe)&  BCi );
            Age =   BCe ^((~BCi)&  BCo );
            Agi =   BCi ^((~BCo)&  BCu );
            Ago =   BCo ^((~BCu)&  BCa );
            Agu =   BCu ^((~BCa)&  BCe );

            Ebe ^= De;
            BCa = ROL(Ebe, 1);
            Egi ^= Di;
            BCe = ROL(Egi, 6);
            Eko ^= Do;
            BCi = ROL(Eko, 25);
            Emu ^= Du;
            BCo = ROL(Emu, 8);
            Esa ^= Da;
            BCu = ROL(Esa, 18);
            Aka =   BCa ^((~BCe)&  BCi );
            Ake =   BCe ^((~BCi)&  BCo );
            Aki =   BCi ^((~BCo)&  BCu );
            Ako =   BCo ^((~BCu)&  BCa );
            Aku =   BCu ^((~BCa)&  BCe );

            Ebu ^= Du;
            BCa = ROL(Ebu, 27);
            Ega ^= Da;
            BCe = ROL(Ega, 36);
            Eke ^= De;
            BCi = ROL(Eke, 10);
            Emi ^= Di;
            BCo = ROL(Emi, 15);
            Eso ^= Do;
            BCu = ROL(Eso, 56);
            Ama =   BCa ^((~BCe)&  BCi );
            Ame =   BCe ^((~BCi)&  BCo );
            Ami =   BCi ^((~BCo)&  BCu );
            Amo =   BCo ^((~BCu)&  BCa );
            Amu =   BCu ^((~BCa)&  BCe );

            Ebi ^= Di;
            BCa = ROL(Ebi, 62);
            Ego ^= Do;
            BCe = ROL(Ego, 55);
            Eku ^= Du;
            BCi = ROL(Eku, 39);
            Ema ^= Da;
            BCo = ROL(Ema, 41);
            Ese ^= De;
            BCu = ROL(Ese, 2);
            Asa =   BCa ^((~BCe)&  BCi );
            Ase =   BCe ^((~BCi)&  BCo );
            Asi =   BCi ^((~BCo)&  BCu );
            Aso =   BCo ^((~BCu)&  BCa );
            Asu =   BCu ^((~BCa)&  BCe );
        }

        //copyToState(state, A)
        state[ 0] = Aba;
        state[ 1] = Abe;
        state[ 2] = Abi;
        state[ 3] = Abo;
        state[ 4] = Abu;
        state[ 5] = Aga;
        state[ 6] = Age;
        state[ 7] = Agi;
        state[ 8] = Ago;
        state[ 9] = Agu;
        state[10] = Aka;
        state[11] = Ake;
        state[12] = Aki;
        state[13] = Ako;
        state[14] = Aku;
        state[15] = Ama;
        state[16] = Ame;
        state[17] = Ami;
        state[18] = Amo;
        state[19] = Amu;
        state[20] = Asa;
        state[21] = Ase;
        state[22] = Asi;
        state[23] = Aso;
        state[24] = Asu;

        #undef    round
}

#include <string.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))


/*************************************************
* Name:        keccak_absorb
* 
* Description: Absorb step of Keccak;
*              non-incremental, starts by zeroeing the state.
*
* Arguments:   - uint64_t *s:             pointer to (uninitialized) output Keccak state
*              - unsigned int r:          rate in bytes (e.g., 168 for SHAKE128)
*              - const unsigned char *m:  pointer to input to be absorbed into s
*              - unsigned long long mlen: length of input in bytes
*              - unsigned char p:         domain-separation byte for different Keccak-derived functions
**************************************************/
static void keccak_absorb(uint64_t *s,
                          unsigned int r,
                          const unsigned char *m, unsigned long long int mlen,
                          unsigned char p)
{
  unsigned long long i;
  unsigned char t[200];

  for (i = 0; i < 25; ++i)
    s[i] = 0;
  
  while (mlen >= r) 
  {
    for (i = 0; i < r / 8; ++i)
      s[i] ^= load64(m + 8 * i);
    
    KeccakF1600_StatePermute(s);
    mlen -= r;
    m += r;
  }

  for (i = 0; i < r; ++i)
    t[i] = 0;
  for (i = 0; i < mlen; ++i)
    t[i] = m[i];
  t[i] = p;
  t[r - 1] |= 128;
  for (i = 0; i < r / 8; ++i)
    s[i] ^= load64(t + 8 * i);
}

/*************************************************
* Name:        keccak_squeezeblocks
* 
* Description: Squeeze step of Keccak. Squeezes full blocks of r bytes each.
*              Modifies the state. Can be called multiple times to keep squeezing,
*              i.e., is incremental.
*
* Arguments:   - unsigned char *h:               pointer to output blocks
*              - unsigned long long int nblocks: number of blocks to be squeezed (written to h)
*              - uint64_t *s:                    pointer to in/output Keccak state
*              - unsigned int r:                 rate in bytes (e.g., 168 for SHAKE128)
**************************************************/
static void keccak_squeezeblocks(unsigned char *h, unsigned long long int nblocks,
                                 uint64_t *s, unsigned int r)
{
  unsigned int i;
  while(nblocks > 0) 
  {
    KeccakF1600_StatePermute(s);
    for(i=0;i<(r>>3);i++)
    {
      store64(h+8*i, s[i]);
    }
    h += r;
    nblocks--;
  }
}

/*************************************************
* Name:        shake128_absorb
* 
* Description: Absorb step of the SHAKE128 XOF.
*              non-incremental, starts by zeroeing the state.
*
* Arguments:   - uint64_t *s:                     pointer to (uninitialized) output Keccak state
*              - const unsigned char *input:      pointer to input to be absorbed into s
*              - unsigned long long inputByteLen: length of input in bytes
**************************************************/
void shake128_absorb(uint64_t *s, const unsigned char *input, unsigned long long inputByteLen)
{
  keccak_absorb(s, SHAKE128_RATE, input, inputByteLen, 0x1F);
}

/*************************************************
* Name:        shake128_squeezeblocks
* 
* Description: Squeeze step of SHAKE128 XOF. Squeezes full blocks of SHAKE128_RATE bytes each.
*              Modifies the state. Can be called multiple times to keep squeezing,
*              i.e., is incremental.
*
* Arguments:   - unsigned char *output:      pointer to output blocks
*              - unsigned long long nblocks: number of blocks to be squeezed (written to output)
*              - uint64_t *s:                pointer to in/output Keccak state
**************************************************/
void shake128_squeezeblocks(unsigned char *output, unsigned long long nblocks, uint64_t *s)
{
  keccak_squeezeblocks(output, nblocks, s, SHAKE128_RATE);
}

/*************************************************
* Name:        shake256
* 
* Description: SHAKE256 XOF with non-incremental API
*
* Arguments:   - unsigned char *output:      pointer to output
*              - unsigned long long outlen:  requested output length in bytes
               - const unsigned char *input: pointer to input
               - unsigned long long inlen:   length of input in bytes
**************************************************/
void shake256(unsigned char *output, unsigned long long outlen, 
              const unsigned char *input,  unsigned long long inlen)
{
  uint64_t s[25];
  unsigned char t[SHAKE256_RATE];
  unsigned long long nblocks = outlen/SHAKE256_RATE;
  size_t i;

  for (i = 0; i < 25; ++i)
    s[i] = 0;
  
  /* Absorb input */
  keccak_absorb(s, SHAKE256_RATE, input, inlen, 0x1F);

  /* Squeeze output */
  keccak_squeezeblocks(output, nblocks, s, SHAKE256_RATE);

  output+=nblocks*SHAKE256_RATE;
  outlen-=nblocks*SHAKE256_RATE;

  if(outlen) 
  {
    keccak_squeezeblocks(t, 1, s, SHAKE256_RATE);
    for(i=0;i<outlen;i++)
      output[i] = t[i];
  }
}


/*************************************************
* Name:        shake256_absorb
*
* Description: Absorb step of the SHAKE256 XOF.
*              non-incremental, starts by zeroeing the state.
*
* Arguments:   - keccak_state *s:   pointer to (uninitialized) output Keccak state
*              - const uint8_t *in: pointer to input to be absorbed into s
*              - size_t inlen:      length of input in bytes
**************************************************/
void shake256_absorb(keccak_state *state, const unsigned char *in, unsigned long long inlen)
{
  keccak_absorb(state->s, SHAKE256_RATE, in, inlen, 0x1F);
}

/*************************************************
* Name:        shake256_squeezeblocks
*
* Description: Squeeze step of SHAKE256 XOF. Squeezes full blocks of
*              SHAKE256_RATE bytes each. Modifies the state. Can be called
*              multiple times to keep squeezing, i.e., is incremental.
*
* Arguments:   - uint8_t *out:    pointer to output blocks
*              - size_t nblocks:  number of blocks to be squeezed
*                                 (written to output)
*              - keccak_State *s: pointer to input/output Keccak state
**************************************************/
void shake256_squeezeblocks(unsigned char *out, unsigned long long nblocks, keccak_state *state)
{
  keccak_squeezeblocks(out, nblocks, state->s, SHAKE256_RATE);
}

/*************************************************
* Name:        keccak_absorb_masked
*
* Description: Absorb step of Keccak;
*              non-incremental, starts by zeroeing the state.
*
* Arguments:   - uint64_t *s: pointer to (uninitialized) output Keccak state
*              - unsigned int r: rate in bytes (e.g., 168 for SHAKE128)
*              - const uint8_t *m: pointer to input to be absorbed into s
*              - size_t mlen: length of input in bytes
*              - uint8_t p: domain-separation byte for different
*                           Keccak-derived functions
**************************************************/
 void keccak_absorb_masked(uint64_t s_masked[25 * (NEWHOPE_MASKING_ORDER + 1)],
                           unsigned int r,
                           const unsigned char* m_masked,
                           unsigned long long mlen,
                           unsigned char p)
{
    size_t i, tmp_mlen;
    unsigned int l;
    unsigned char t[200*(NEWHOPE_MASKING_ORDER+1)] = { 0 };
    tmp_mlen = mlen;
    /* Zero state */
    for (l = 0; l < (NEWHOPE_MASKING_ORDER + 1); l++)
    {
        for (i = 0; i < 25; i++)
            s_masked[i+(25*l)] = 0;
    }
    
    while (mlen >= r) {
        for (l = 0; l < (NEWHOPE_MASKING_ORDER + 1); l++)
        {
            for (i = 0; i < r / 8; i++)
                s_masked[i + (l * 25)] ^= load64(m_masked + (8 * i) + (l*tmp_mlen));
        }
        KeccakF1600_StatePermute_masked(s_masked);
        mlen -= r;
        m_masked += r;
    }
    for (l = 0; l < (NEWHOPE_MASKING_ORDER + 1); l++)
    {
        for (i = 0; i < mlen; i++)
            t[i + (200 * l)] = m_masked[i + (l * tmp_mlen)];
        t[i] = p; //Only in the first share
        t[r - 1] |= 128; //Only in the first share
        for (i = 0; i < r / 8; i++)
            s_masked[i + (25 * l)] ^= load64(t + (8 * i) + (200 * l));
    }
}

/*************************************************
* Name:        keccak_squeezeblocks
*
* Description: Squeeze step of Keccak. Squeezes full blocks of r bytes each.
*              Modifies the state. Can be called multiple times to keep
*              squeezing, i.e., is incremental.
*
* Arguments:   - uint8_t *h: pointer to output blocks
*              - size_t nblocks: number of blocks to be squeezed (written to h)
*              - uint64_t *s: pointer to input/output Keccak state
*              - unsigned int r: rate in bytes (e.g., 168 for SHAKE128)
**************************************************/
static void keccak_squeezeblocks_masked(unsigned char* out_masked,
                                        size_t nblocks,
                                        uint64_t s_masked[25*(NEWHOPE_MASKING_ORDER+1)],
                                        unsigned int r,
                                        size_t outlen)
{
    unsigned int i,l,offset;
    if (outlen == 0)
        offset = r * nblocks;
    else
        offset = outlen;
    while (nblocks > 0)
    {
        KeccakF1600_StatePermute_masked(s_masked);
        for (l = 0; l < (NEWHOPE_MASKING_ORDER + 1); l++)
        {
            for (i = 0; i < r / 8; i++)
                store64(out_masked + (8 * i) + (l* offset), s_masked[i+(l*25)]);
        }
        out_masked += r;
        --nblocks;
    }
    return;
}


/*************************************************
* Name:        shake256_absorb
*
* Description: Absorb step of the SHAKE256 XOF.
*              non-incremental, starts by zeroeing the state.
*
* Arguments:   - keccak_state *s:   pointer to (uninitialized) output Keccak state
*              - const uint8_t *in: pointer to input to be absorbed into s
*              - size_t inlen:      length of input in bytes
**************************************************/
void shake256_absorb_masked(keccak_state_masked* state_masked, const unsigned char* in_masked, size_t inlen)
{
    keccak_absorb_masked(state_masked->s_masked, SHAKE256_RATE, in_masked, inlen, 0x1F);
}

/*************************************************
* Name:        shake256_squeezeblocks
*
* Description: Squeeze step of SHAKE256 XOF. Squeezes full blocks of
*              SHAKE256_RATE bytes each. Modifies the state. Can be called
*              multiple times to keep squeezing, i.e., is incremental.
*
* Arguments:   - uint8_t *out:    pointer to output blocks
*              - size_t nblocks:  number of blocks to be squeezed
*                                 (written to output)
*              - keccak_State *s: pointer to input/output Keccak state
**************************************************/
void shake256_squeezeblocks_masked(unsigned char* out_masked, size_t nblocks, keccak_state_masked* state_masked, size_t outlen)
{
    keccak_squeezeblocks_masked(out_masked, nblocks, state_masked->s_masked, SHAKE256_RATE, outlen);
}

void shake256_masked(unsigned char* out_masked, size_t outlen, const unsigned char* in_masked, size_t inlen)
{
    unsigned int i,l;
    size_t tmp_outlen;
    size_t nblocks = outlen / SHAKE256_RATE;
    unsigned char t_masked[SHAKE256_RATE*(NEWHOPE_MASKING_ORDER+1)];
    keccak_state_masked state_masked;

    tmp_outlen = outlen;
    shake256_absorb_masked(&state_masked, in_masked, inlen);
    shake256_squeezeblocks_masked(out_masked, nblocks, &state_masked, tmp_outlen);

    out_masked += nblocks * SHAKE256_RATE;
    outlen -= nblocks * SHAKE256_RATE;

    if (outlen) {
        shake256_squeezeblocks_masked(t_masked, 1, &state_masked, 0);
        for (l = 0; l < (NEWHOPE_MASKING_ORDER + 1); l++)
        {
            for (i = 0; i < outlen; i++)
                out_masked[i + (l * tmp_outlen)] = t_masked[i + (l * SHAKE256_RATE)];
        }
    }
}