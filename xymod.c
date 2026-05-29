#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xymod.h"

// Neighbour lookup table

void ntable(const int L, const int nbr[L*L][NBS])
{
	/************* HACK ALERT!! We un-const nbr → ncnbr *************/
	typedef int (*tbl)[L][NBS];
	tbl ncnbr = (tbl)nbr;
	/****************************************************************/
	const int L1 = L-1;
	for (int x = 0; x < L; ++x) {
		for (int y = 0; y < L; ++y) {
			ncnbr[x][y][0] = L*x           + ((y+1)%L);   // right
			ncnbr[x][y][1] = L*((x+1)%L)   + y;           // down
			ncnbr[x][y][2] = L*x           + ((y+L1)%L);  // left
			ncnbr[x][y][3] = L*((x+L1)%L ) + y;           // up
		}
	}
}

// A single Langevin (Euler-Maruyama) lattice update

inline void Langevin_update(const int N, double hnew[N], const double h[N], const int nbr[N][NBS], const double sig, const double Jdt, rng_t* const prng)
{
	// fill new spin buffer with Wiener noise
	spin_normal(N,hnew,sig,prng);

	for (int i = 0; i<N; ++i) {

		// calculate torque
		const double hi   = h[i];
		const int* const nbri = nbr[i];
		double torqi = sin(hi - h[nbri[0]]) + sin(hi - h[nbri[1]]) + sin(hi - h[nbri[2]]) + sin(hi - h[nbri[3]]);

		// Euler-Maruyama update
		hnew[i] += hi - Jdt*torqi;
	}
}

// A single Wolff lattice update

inline int Wolff_update(int N, uvec_t v[N], const int nbr[N][NBS], bool cluster[N], int stack[N], double Jbeta, rng_t* const prng)
{
	// Random reflection axis
	const double angle = 2.0*M_PI*rng_rand(prng);
	const double rx = cos(angle);
	const double ry = sin(angle);

	memset(cluster,0,(size_t)N*sizeof(bool));
//	for (int i = 0; i < N; ++i) cluster[i] = false;

	int i0 = (int)((double)N*rng_rand(prng));
	stack[0]    = i0;
	int top     = 1;
	cluster[i0] = true;

	// Build cluster
	while (top > 0) {

		--top;
		const int i = stack[top];

		const double sproj = v[i].x*rx + v[i].y*ry;

		for (int k = 0; k < NBS; ++k) { // for each neighbour

			const int in = nbr[i][k];

			if (cluster[in]) continue; // already there

			const double tproj = v[in].x*rx + v[in].y*ry;

			if (sproj*tproj > 0.0) {

				const double pij = 1.0 - exp(-2.0*Jbeta*sproj*tproj);
				if (rng_rand(prng) < pij) { // add to cluster
					cluster[in]  = true;
					stack[top++] = in;
				}
			}
		}
	}

    // Reflect cluster
    int c = 0;
	for (int i = 0; i < N; ++i) {
		if (cluster[i]) {
			const double dot = v[i].x*rx + v[i].y*ry;
			v[i].x -= 2.0*dot*rx;
			v[i].y -= 2.0*dot*ry;
			++c;
		}
	}

	return c;
}

inline int quadrant(const double x, const double y)
{
    if (x > 0.0) return y >= 0.0 ? 0 : 3;
    if (x < 0.0) return y >  0.0 ? 1 : 2;
    if (y > 0.0) return 1; // x = 0, y > 0
    return 3;              // x = 0, y < 0
}

vlist_t* vortices(const int N, const uvec_t v[N], const int nbr[N][NBS], int* const vcount)
{
	// IMPORTANT: to avoid memory leaks, after use you MUST run vlist_free() on the
	// vlist_t pointer returned by this function to free the vortex list memory!!!

	*vcount = 0;
	vlist_t* vlist = NULL;
	for (int i = 0; i < N; ++i) {
		// around plaque
		const int j = nbr[i][0]; // right
		const int k = nbr[j][1]; // right down
		const int l = nbr[i][1]; // down
		double X[4], Y[4];
		X[0] = v[i].x; Y[0] = v[i].y;
		X[1] = v[j].x; Y[1] = v[j].y;
		X[2] = v[k].x; Y[2] = v[k].y;
		X[3] = v[l].x; Y[3] = v[l].y;
		int dqsum = 0;
		int qprev = quadrant(X[3],Y[3]);
		for (int q = 0; q < 4; ++q) {
			const int qcurr = quadrant(X[q],Y[q]);
			int dq = qcurr-qprev;
			// wrap-around
			if      (dq == +3) dq = -1;
			else if (dq == -3) dq = +1;
			else if (dq == +2 || dq == -2) { // resolve diagonal jumps (Axiom 1 violation)
				// Use the cross product of S_{q-1} and S_i S_prev is at index (q-1+n)%n
				const int iprev = q == 0 ? 3 : q - 1; // wrap
				dq = X[iprev]*Y[q] > Y[iprev]*X[q] ? 2 : -2;
			}
			dqsum += dq;
			qprev = qcurr;
		}
		const int w = dqsum/4;
		if (w !=  0){
			 ++(*vcount);
			vlist = vlist_push(vlist,w,i);
		 }
	}
	return vlist;
}

void phase_grad(const int N, const uvec_t v[N], const int nbr[N][NBS], vec_t dv[N])
{
	for (int i = 0; i < N; ++i) {
		const uvec_t vi = v[i];
		const uvec_t vx = v[nbr[i][0]]; // right
		const uvec_t vy = v[nbr[i][1]]; // down
		dv[i].x =  vi.x*(vx.y - vi.y) - vi.y*(vx.x - vi.x);
		dv[i].y =  vi.x*(vy.y - vi.y) - vi.y*(vy.x - vi.x);
	}
}
