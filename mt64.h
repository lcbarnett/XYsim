#ifndef MT64_H
#define MT64_H

// 64-bit Mersenne Twister PRNG (thread-safe)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <sys/random.h>

#define NN 312
#define MM 156
#define MATRIX_A UINT64_C(0xB5026F5AA96619E9)
#define UM UINT64_C(0xFFFFFFFF80000000) // Most  significant 33 bits
#define LM UINT64_C(0x7FFFFFFF)         // Least significant 31 bits

typedef uint64_t ruint_t;

typedef struct {
	ruint_t mt[NN];
	int      mti;
	int      iset;
	double   gset;
} rng_t;

// SplitMix64: generate next 64-bit seed (make sure to set 'state' to the
// master seed first). This is a half-decent PRNG in its own right ;-)
static inline ruint_t rng_sm64(ruint_t* const state)
{
	ruint_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
	z = (z ^ (z >> 30)) *  UINT64_C(0xBF58476D1CE4E5B9);
	z = (z ^ (z >> 27)) *  UINT64_C(0x94D049BB133111EB);
	return z ^ (z >> 31);
}

// get a seed from /dev/urandom
static inline ruint_t rng_getrandom()
{
	ruint_t seed;
	// this call shouldn't fail :-O
	if (getrandom(&seed,sizeof(ruint_t),0) != sizeof(ruint_t)) fprintf(stderr,"WARNING: getrandom() failed\n");
	return seed;
}

// initializes state
static inline ruint_t rng_seed(rng_t* const pstate, ruint_t seed)
{
	// The generator MUST be seeded before use!!! If 0 is supplied we read /dev/urandom.

	if (seed == 0) seed = rng_getrandom();
	const ruint_t eseed = seed; // effective seed (to be returned for reproducability)

	pstate->mt[0] = rng_sm64(&seed); // we don't trust "small" user-supplied seeds, so we munge all seeds through SplitMix64 :-)
	for (pstate->mti=1; pstate->mti<NN; pstate->mti++) pstate->mt[pstate->mti] = (UINT64_C(6364136223846793005) * (pstate->mt[pstate->mti-1] ^ (pstate->mt[pstate->mti-1] >> 62)) + (ruint_t)pstate->mti);

	pstate->iset = 0;
	pstate->gset = 0.0;

    return eseed; // you may use the returned seed for reproducability (even if 0 was previously supplied)
}

// copy full state of rng (for save/restore)
static inline void rng_copy(rng_t* const dstate, const rng_t* const sstate)
{
	memcpy(dstate->mt,sstate->mt,NN*sizeof(ruint_t));
	dstate->mti  = sstate->mti;
	dstate->iset = sstate->iset;
	dstate->gset = sstate->gset;
}

// generates a random number on [0, 2^64-1]-interval
static inline ruint_t rng_uint(rng_t* const pstate)
{
    int i;
    ruint_t x;
    static const ruint_t mag01[2]={UINT64_C(0), MATRIX_A};

    if (pstate->mti >= NN) { // generate NN words at one time

        // if init_genrand64() has not been called, a default initial seed is used
        // if (pstate->mti == NN+1) mt64_seed(pstate,UINT64_C(5489));
        // No! We don't want this. The user MUST seed first! (probably segfault if they forget)

        for (i=0;i<NN-MM;i++) {
            x = (pstate->mt[i]&UM)|(pstate->mt[i+1]&LM);
            pstate->mt[i] = pstate->mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        for (;i<NN-1;i++) {
            x = (pstate->mt[i]&UM)|(pstate->mt[i+1]&LM);
            pstate->mt[i] = pstate->mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        x = (pstate->mt[NN-1]&UM)|(pstate->mt[0]&LM);
        pstate->mt[NN-1] = pstate->mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];

        pstate->mti = 0;
    }

    x = pstate->mt[pstate->mti++];

    x ^= (x >> 29) & UINT64_C(0x5555555555555555);
    x ^= (x << 17) & UINT64_C(0x71D67FFFEDA60000);
    x ^= (x << 37) & UINT64_C(0xFFF7EEE000000000);
    x ^= (x >> 43);

    return x;
}

// NOTE: IEEE 754 double-precision gives 53-bit precision (2^53 = 9007199254740992).
// The high bits of the MT64 are considered "more random", so for double-precision fp
// we right-shift by 11 bits, to use just the 53 high bits :-)

// generates a random number uniform on [0,1)-real-interval
static inline double rng_rand(rng_t* const pstate)
{
	return (double)(rng_uint(pstate) >> 11) * (1.0/9007199254740992.0);
}

// generates a random number uniform on (0,1)-real-interval
static inline double rng_rand_pos(rng_t* const pstate)
{
	double x;
	do x = (double)(rng_uint(pstate) >> 11) * (1.0/9007199254740992.0); while (x == 0.0);
	return x;
}

// return integer in range 0,...,range-1 (be careful of overflow!!)
#define RANDI(itype,range,pstate) (itype)((double)(range)*rng_rand(pstate));

// generates a normally distributed random number from N(0,1)
static inline double rng_randn(rng_t* const pstate)
{
    if (pstate->iset) {
	    pstate->iset=0;
	    return pstate->gset;
    }

    double v1,v2,rsq;
    do {
	    v1 = 2.0*rng_rand(pstate)-1.0;
	    v2 = 2.0*rng_rand(pstate)-1.0;
	    rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    const double fac = sqrt(-2.0*log(rsq)/rsq);
    pstate->gset = fac*v1;
    pstate->iset = 1;
    return fac*v2;
}

#undef NN
#undef MM
#undef MATRIX_A
#undef UM
#undef LM

#endif // MT64_H
