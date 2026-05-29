#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>

#include "mt64.h"
#include "vlist.h"

// Number of lattice neighbours

#define NBS (4)

// Berezinskii-Kosterlitz-Thouless phase transition temperature

static const double TBKT = 0.89290;

// 2D vector struct type

typedef struct {
	double x;
	double y;
} vec_t;

typedef vec_t uvec_t; // unit vector is the same type, but we use a different name for clairty

// Build neighbour lookup table

void ntable(const int L, const int nbr[L*L][NBS]);

// Swap double pointers (for buffer-swap trick) – usual caveats for UB

#define dpswap(x1,x2) {double* const xtmp = x1; x1 = x2; x2 = xtmp;}

// LATTICE UPDATES /////////////////////////////////////////////////////////////

// A single Langevin (Euler-Maruyama) lattice update

void Langevin_update(const int N, double hnew[N], const double h[N], const int nbr[N][NBS], const double sig, const double Jdt, rng_t* const prng);

// A single Wolff lattice update

int Wolff_update(int N, uvec_t v[N], const int nbr[N][NBS], bool cluster[N], int stack[N], double Jbeta, rng_t* const prng);

// Multiple Wolff lattice updates to achieve a minumum number of spin reflections

static inline int Wolff_updatex(int N, uvec_t v[N], const int nbr[N][NBS], bool cluster[N], int stack[N], double Jbeta, const int minflips, rng_t* const prng)
{
	// Enough updates that at least minflips sites are reflected
	int updates = 0;
	for (int flips = 0; flips < minflips; ++updates) flips += Wolff_update(N,v,nbr,cluster,stack,Jbeta,prng);
	return updates;
}

// Find vortices

vlist_t* vortices   (const int N, const uvec_t v[N], const int nbr[N][NBS], int* const vcount);

// Calculate phase gradient

void     phase_grad (const int N, const uvec_t v[N], const int nbr[N][NBS], vec_t dv[N]);

// Utilities

static inline void angle2uvec(const int n, uvec_t v[n], const double h[n])
{
	for (int i = 0; i < n; ++i) {
		v[i].x = cos(h[i]);
		v[i].y = sin(h[i]);
	}
}

static inline void uvec2angle(const int n, double h[n], const uvec_t v[n])
{
	for (int i = 0; i < n; ++i) h[i] = atan2(v[i].y,v[i].x);
}

static inline void uvec_uniform(const int n, uvec_t v[n], rng_t* const prng)
{
	// Fill unit vector buffer iid from a uniform distribution of angle on [0,2π).
	for (int i = 0; i < n; ++i) {
		const double h = 2.0*M_PI*rng_rand(prng);
		v[i].x = cos(h);
		v[i].y = sin(h);
	}
}

static inline void spin_normal(const int n, double h[n], const double sig, rng_t* const prng)
{
	// Fill spin buffer iid from a normal distribution N(0,sig^2).
	for (int i = 0; i < n; ++i) h[i] = sig*rng_randn(prng);
}
