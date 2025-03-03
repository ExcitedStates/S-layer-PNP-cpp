#include <iostream>
#include <fstream>
#include <sstream>
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
	int n_save = 0;
	ostringstream number_stream;

	/* check args */
	if (argc < 2) {
		cout << "No input file!" << endl;
		return 1;
	} else {
		if (argc < 3) {
			cout << "No time step number provided, will only run 1 step." << endl;
		} else {
			n_t = atoi( argv[2] );
			if (argc > 3) {
				n_save = atoi( argv[3] );
			} else {
				n_save = 100;
			}
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
	double *km;

	double *bulk;
	double *amo_z, *amo_r;
	double *slp;

	double *n_display;
	double *t_after;

	int m0, n0, l0, m1, n1, l1, m2, n2, l2, m3, n3, l3;
	int t_offset;

	double *CNew;
	double *EzNew, *ErNew, *EqNew;
	double *JzNew, *JrNew, *JqNew;

	time_t t0, t1, t2;

	/* parse profile dimensions */
	string filename = argv[1];
	readBinHeader( filename, &t_offset,
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
	km = new double;
	bulk = new double[ 1 * n0 * l0 * n_i ];
	amo_z = new double;
	amo_r = new double;
	slp = new double[3];
	n_display = new double;
	t_after = new double;

	/* read all variables */
	readBin( argv[1], n_i, C, Ez, Er, Eq, Jz, Jr, Jq,
		&m0, &n0, &l0, &m1, &n1, &l1, &m2, &n2, &l2, &m3, &n3, &l3, z, d_m,
	 	dqq, dx, dt, R, a, b, c, km, bulk , amo_z, amo_r, slp, n_display, t_after );

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
			dqq, dx, dt, R, a, b, c, km, bulk, amo_z, amo_r, slp, n_display, t_after
		 	);
		*t_after -= 1.0 * *n_display;
		t2 = std::time(0);
		t_eta = 0.1 * round( 10 * (t2-t0)/(t+1.0)*(n_tt-t-1.0));
		cout << "step " << ((t+1) * *n_display) << " done in " << 0.1*round(10*(t2-t1)) << " s. ";
		cout << t_eta << "s to go..." << endl;

		/* save intermediate result */
		if ((t+1) % n_save == 0) {
			if ((t+1) == n_tt) {
				cout << "saving final result... ";
			} else {
				cout << "saving intermediate result... ";
			}
			t_offset += n_save * *n_display;
			/* write all variables */
			number_stream.str("");
			number_stream << t_offset;
			filename = filename.substr( 0, filename.length()-(5+(int)floor(log10(t_offset))) ) + number_stream.str() + ".pnl";
			writeBin( filename, n_i, t_offset,
				C, Ez, Er, Eq, Jz, Jr, Jq,
				&m0, &n0, &l0, &m1, &n1, &l1, &m2, &n2, &l2, &m3, &n3, &l3 );
			cout << "done!" << endl;
		}
	}

	cout << "simulation concluded" << endl;

	return 0;

}
