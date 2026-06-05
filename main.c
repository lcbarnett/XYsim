#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <raylib.h>

#include "xymod.h"
#include "clap.h"

#define WTEXTMAX 200

// See 'fonts' directory. Currently contains
//
// DejaVuSansMono.ttf
// LiberationMono-Regular.ttf
// NotoMono-Regular.ttf

// Raylib verbosity flags - see SetTraceLogLevel()
//
// LOG_ALL:     0
// LOG_TRACE:   1
// LOG_DEBUG:   2
// LOG_INFO:    3
// LOG_WARNING: 4
// LOG_ERROR:   5
// LOG_FATAL:   6
// LOG_NONE:    7

static inline void help_message(const double TBKT)
{
	printf("\n\033[4mInteraction\033[0m\n"); // underline!
	printf("\n*  The 'A' key toggles the update algorithm between Langevin and Wolff.\n");
	printf("\n*  The 'G' key toggles animation between spin field and spin-gradient field.\n");
	printf("\n*  The 'C' key sets the temperature to the BKT critical value T ≈ %.4f.\n",TBKT);
	printf("\n*  The 'I' key re-initialises the lattice state to uniform random spins.\n");
	printf("\n*  RIGHT/LEFT arrow keys increase/decrease the updates-per-frame rate.\n");
	printf("\n*  UP/DOWN arrow keys increase/decrease temperature by 'Tinc' (default 0.05)\n");
	printf("\n*  The SPACE key pauses/resumes the simulation.\n");
	printf("\n*  The ESC key exits the program.\n\n");
}

static inline void set_window_text(char* const wtext, const bool lgv, const bool grad, const double T, const int upf, const size_t itr, const int vnum)
{
	snprintf(wtext,WTEXTMAX,"%s field : T = %.4f : %s x %d : updates %zu : vortices = %d",grad?"Spin-gradient":"Spin",T,lgv?"Langevin":"Wolff",upf,itr,vnum);
}

int main(int argc, char* argv[])
{
	--argc; ++argv; // strip program name (for CLAP_* macros)

	if (argc == 1 && strcasecmp(argv[0],"-help") == 0) {
		help_message(TBKT);
		return EXIT_SUCCESS;
	}

	printf("\n\033[4mParameters\033[0m (override defaults as command-line switches)\n\n"); // underline!
	CLAP_VARG(lgv,        bool,     true,             "Langevin update (else Wolff)");
	CLAP_VARG(upf,        int,      1,                "Langevin/Wolff updates per frame");
	CLAP_VARG(grad,       bool,     false,            "Display spin gradient field rather than spin field");
	CLAP_CARG(L,          int,      32,               "Lattice length");
	CLAP_CARG(FBI,        int,      0,                "Burn-in flips (x N)");
	CLAP_VARG(T,          double,   TBKT,             "Temperature (default: BKT critical value)");
	CLAP_CARG(Tinc,       double,   0.05,             "Temperature increment for interactive adjustment");
	CLAP_CARG(dt,         double,   0.01,             "Langevin Euler-Maruyama integration step size");
	CLAP_CARG(scr,        int,      1364,             "Raylib window size (pixels)");
	CLAP_CARG(fps,        int,      60,               "Raylib target (maximum) frames per second");
	CLAP_CARG(alength,    float,    0.8f,             "Raylib quiver plot arrow shaft length");
	CLAP_CARG(athick,     float,    1.0f,             "Raylib quiver plot arrow shaft thickness");
	CLAP_CARG(aheads,     float,    0.2f,             "Raylib quiver plot arrow head size");
	CLAP_CARG(hangle,     float,    0.5f,             "Raylib quiver plot arrow head angle (radians)");
	CLAP_CARG(carrow,     bool,     true,             "Raylib quiver plot arrow - center at lattice site");
	CLAP_CARG(vcrad,      float,    3.5f,             "Raylib vortex circle radius (pixels)");
	CLAP_CARG(fsize,      int,      20,               "Raylib font size");
	CLAP_CARG(rseed,      ulong,    0,                "PRNG seed (0 for random random seed)");

	help_message(TBKT);

	rng_t rng;
	const ruint_t aseed = rng_seed(&rng,rseed);
	printf("Actual seed = %lu\n\n",aseed);

	const int N = L*L;

	double beta = 1.0/T;
	double sig  = sqrt(2.0*T*dt); // Wiener noise intensity

	// Neighbour table
	const int nbr[N][NBS];
	ntable(L,nbr);

    // Wolff update cluster and stack
	bool cluster[N];
	int  stack[N];

	// Langevin update spin angle buffers
	double h1[N];
	double h2[N];

	// Spin angle buffer pointers (for pointer-swap trick)
	double* h    = h1;
	double* hold = h2;

	// Spin (unit) vector buffer
	uvec_t v[N];

	// Gradient vector buffer
	vec_t dv[N];

	// Spin or spin-gradient
	const vec_t* vf = (grad ? dv : v);

	// Raylib configuration (screen metrics, colours, initialisation, etc.)

	const int ippc = scr/(L+1); // pixels per cell x

	const int   iyoff = 20; // extra vertical pixels for text
	const float yoff  = (float)iyoff;

	const int wx = (L+1)*ippc;         // window x (pixels)
	const int wy = (L+1)*ippc + iyoff; // window y (pixels)

	const float ppc = (float)ippc;     // pixels per cell x

	const float asize = alength*ppc;   // scaled arrow length
	const float ahead = aheads*ppc;    // scaled arrow head size

	float latx[N], laty[N]; // lattice coordinates
	for (int i = 0; i < N; ++i) latx[i] = ppc*((float)(i/L)+1.0f);
	for (int i = 0; i < N; ++i) laty[i] = ppc*((float)(i%L)+1.0f) + yoff;

	const Color lcol = ColorFromHSV(360.0f, 0.0f, 0.0f); // arrow line colour
	const Color vcol = ColorFromHSV(  0.0f, 1.0f, 1.0f); // vortex colour
	const Color acol = ColorFromHSV(255.0f, 1.0f, 1.0f); // anti-vortex colour
	const Color tcol = ColorFromHSV(127.0f, 0.5f, 0.0f); // text colour

	const float vrad = (L < 47 ? 1.25f*vcrad : vcrad);   // bit of a hack :-O

	char wtext[WTEXTMAX+1];                              // window text buffer

	// Initialise Raylib window
	SetConfigFlags(FLAG_MSAA_4X_HINT);
	SetTraceLogLevel(LOG_WARNING); // set verbosity level
 	InitWindow(wx,wy,"XY model");
	SetTargetFPS(fps);

	// Load font (after initialising window!) - see top of this file
	const char fontf[] = "fonts/mono_noto.ttf";
	const Font font = LoadFontEx(fontf,fsize,NULL,0);
	if (font.texture.id == 0) TraceLog(LOG_WARNING,"Failed to load font! Falling back to default.");

	// Main animation loop

    size_t itr    = 0;
    bool   paused = false;
	int    vnum   = 0;

    while (!WindowShouldClose()) { // ESC key exits loop

		BeginDrawing();

		// Update spin lattice

		if (!paused) {

			if (itr == 0) {
				// initialise lattice to uniform distribution on [0,2π) and run burn-in updates
				uvec_uniform(N,v,&rng);
				Wolff_updatex(N,v,nbr,cluster,stack,beta,N*FBI,&rng); // burn-in
				if (lgv) uvec2angle(N,h,v); // spin angle form needed for subsequent Langevin_update()
			}
			else {
				// update lattice (Langevin or Wolff)
				if (lgv) {
					for (int u = 0; u < upf; ++u) {
						dpswap(h,hold);
						Langevin_update(N,h,hold,nbr,sig,dt,&rng);
					}
					angle2uvec(N,v,h); // vector form needed for phase_grad() and vcalc_vortices()
				}
				else {
					for (int u = 0; u < upf; ++u) Wolff_update(N,v,nbr,cluster,stack,beta,&rng);
				}
			}

			// calculate phase gradient if requested
			if (grad) phase_grad(N,v,nbr,dv);

			// calculate vortices
			vlist_t* const vlist = vortices(N,v,nbr,&vnum);

			ClearBackground(WHITE);

			// draw window text
			set_window_text(wtext,lgv,grad,T,upf,itr,vnum);
			DrawTextEx(font,wtext,(Vector2){5.0f,3.0f},(float)fsize,1.0f,tcol);

			// draw arrows
			for (int i = 0; i < N; ++i) {

				// lattice site coordinates (pixels)
				const float x = latx[i];
				const float y = laty[i];

				// arrow (scaled spin/gradient vector; pixels)
				const float vx = asize*(float)vf[i].x;
				const float vy = asize*(float)vf[i].y;

				// arrow origin (pixels) - centre on lattice site if requested
				const float ox = (carrow ? x - vx/2.0f : x);
				const float oy = (carrow ? y - vy/2.0f : y);

				// arrow terminus (pixels)
				const float tx = ox + vx;
				const float ty = oy + vy;

				// spin vector angle
				const float angle = atan2f(vy,vx);

				// arrowhead vectors (triangle; third vertex is terminus of spin vector)
				const float ah1x = tx-ahead*cosf(angle+hangle);
				const float ah1y = ty-ahead*sinf(angle+hangle);
				const float ah2x = tx-ahead*cosf(angle-hangle);
				const float ah2y = ty-ahead*sinf(angle-hangle);

				// draw the arrow shaft
				DrawLineEx((Vector2){ox,oy},(Vector2){tx,ty},athick,lcol);

				// draw the arrow head
				DrawTriangle((Vector2){tx,ty},(Vector2){ah1x,ah1y},(Vector2){ah2x,ah2y},lcol);
			}

			// draw vortices
			for (const vlist_t* p = vlist; p; p = p->next) { // iterate through vortex list

				// vortex spec
				const int i = p->vidx; // lattice site index
				const int w = p->wnum; // winding number (±1)

				// circle centre (pixels)  +0.5*ppc  to put in centre of plaque
				const float x = latx[i] + 0.5f*ppc;
				const float y = laty[i] + 0.5f*ppc;
				DrawCircleV((Vector2){x,y},vrad,w > 0 ? acol : vcol);
			}

			// free vortex list
			vlist_free(vlist);

			++itr;
		}

		// process keypresses
		int key;
		while ((key = GetKeyPressed()) != 0) { // note keypresses are queued up
			switch (key) {
				case KEY_A: // toggle update between Langevin and Wolff
					lgv = !lgv;
					if (lgv) uvec2angle(N,h,v); // spins for Langevin (already got vectors for Wolff)
					break;
				case KEY_G: // toggle spin/spin-gradient field
					grad = !grad;
					vf = (grad ? dv : v); // the vector field to animate (spin or spin-gradient
					break;
				case KEY_I: // re-initialise lattice to uniform random
					uvec_uniform(N,v,&rng);
					if (lgv) uvec2angle(N,h,v); // spin angle form needed for subsequent Langevin_update()
					break;
				case KEY_SPACE: // pause animation
					paused = !paused;
					break;
				case KEY_LEFT: // decrease updates per iteration (unless at zero)
					if (upf > 1) --upf; else fprintf(stderr,"*** Cannot reduce updates-per-frame further!\n");
					break;
				case KEY_RIGHT: // increase updates per iteration
					++upf;
					break;
				case KEY_UP: // increase temperature
					T   += Tinc;
					beta = 1.0/T;
					sig  = sqrt(2.0*T*dt);
					break;
				case KEY_DOWN: // decrease temperature
					if (T >= Tinc) {
						T   -= Tinc;
						beta = 1.0/T;
						sig  = sqrt(2.0*T*dt);
					}
					break;
				case KEY_C: // BKT critical temperature
					T    = TBKT;
					beta = 1.0/T;
					sig  = sqrt(2.0*T*dt);
					break;
				default:
					fprintf(stderr,"*** Unhandled key\n");
			}
			set_window_text(wtext,lgv,grad,T,upf,itr,vnum);
			DrawRectangleV((Vector2){0.0f,0.0f},(Vector2){1000.0f,25.0f},WHITE); // clear text drawing space
			DrawTextEx(font,wtext,(Vector2){5.0f,3.0f},(float)fsize,1.0f,tcol);
		}

		if (!paused) {
			snprintf(wtext,WTEXTMAX,"FPS = %3d",GetFPS());
			DrawTextEx(font,wtext,(Vector2){(float)wx-110.0f,3.0f},(float)fsize,1.0f,tcol);
		}

		EndDrawing();
	}

	// Clean up

	UnloadFont(font);
    CloseWindow();

	return EXIT_SUCCESS;
}
