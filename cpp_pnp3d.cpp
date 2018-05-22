#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>
#include <stdlib.h>

using namespace std;

#include "binaryIO.cpp"
#include "solver_pnp3d.cpp"

const int n_i = 3;

int main(int argc, char *argv[]) {

	int n_t = 1;

	/* check args */
	if (argc < 2) {
		cout << "No input file!" << endl;
		return 1;
	} else {
		if (argc < 3) {
			cout << "No time step number provided, will only run 1 step." << endl;
		} else {
			n_t = atoi( argv[2] );
		}
	}

	/* variable declarations */
	double *C;
	double *Ez, *Er, *Eq;
	double *Jz, *Jr, *Jq;

	double *z, *d_m;

	double *dqq, *dx, *dt;
	double *R;
	double *a, *b, *c;

	double *bulk;
	double *amo_z, *amo_r;
	double *slp;

	double *n_display;

	int m0, n0, l0, m1, n1, l1, m2, n2, l2, m3, n3, l3;
	int t_offset;

	double *CNew;
	double *EzNew, *ErNew, *EqNew;
	double *JzNew, *JrNew, *JqNew;

	time_t t0, t1, t2;

	/* parse profile dimensions */
	readBinHeader( argv[1], &t_offset,
		&m0, &n0, &l0, &m1, &n1, &l1, &m2, &n2, &l2, &m3, &n3, &l3 );

	/* declare arrays */
	C = new double[ m0 * n0 * l0 * n_i ];
	Ez = new double[ m1 * n1 * l1 ];
	Er = new double[ m2 * n2 * l2 ];
	Eq = new double[ m3 * n3 * l3 ];
	Jz = new double[ m1 * n1 * l1 * n_i ];
	Jr = new double[ m2 * n2 * l2 * n_i ];
	Jq = new double[ m3 * n3 * l3 * n_i ];
	z = new double[ n_i ];
	d_m = new double[ n_i ];
	dqq = new double;
	dx = new double;
	dt = new double;
	R = new double[ m0 * n0 * (l0+1) ];
	a = new double;
	b = new double;
	c = new double;
	bulk = new double[ 1 * n0 * l0 * n_i ];
	amo_z = new double;
	amo_r = new double;
	slp = new double[3];
	n_display = new double;

	/* read all variables */
	readBin( argv[1], n_i, C, Ez, Er, Eq, Jz, Jr, Jq,
		&m0, &n0, &l0, &m1, &n1, &l1, &m2, &n2, &l2, &m3, &n3, &l3, z, d_m,
	 	dqq, dx, dt, R, a, b, c, bulk, amo_z, amo_r, slp, n_display );

	/* do the PNP business */
	if (n_t < *n_display) {
		n_display[0] = 1.0 * n_t;
	}
	int n_tt = ceil( n_t / *n_display );

	double t_eta = 0;
	t0 = std::time(0);
	for (int t = 0; t < n_tt; t++) {
		t1 = std::time(0);
		solver_pnp3d(C, Ez, Er, Eq, Jz, Jr, Jq,
			m0, n0, l0, m1, n1, l1, m2, n2, l2, m3, n3, l3, z, d_m,
			dqq, dx, dt, R, a, b, c, bulk, amo_z, amo_r, slp, n_display
		 	);
		t2 = std::time(0);
		t_eta = 0.1 * round( 10 * (t2-t0)/(t+1.0)*(n_tt-t-1.0));
		cout << "step " << ((t+1) * *n_display) << " done in " << 0.1*round(10*(t2-t1)) << " s. ";
		cout << t_eta << "s to go..." << endl;
	}

	t_offset += n_t;
	/* write all variables */
	string filename = "results.pnl";
	writeBin( filename, n_i, t_offset,
		C, Ez, Er, Eq, Jz, Jr, Jq,
		&m0, &n0, &l0, &m1, &n1, &l1, &m2, &n2, &l2, &m3, &n3, &l3 );

	return 0;

}
