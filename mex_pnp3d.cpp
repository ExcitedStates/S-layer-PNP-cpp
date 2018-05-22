#include <mex.h>
#include <matrix.h>
#include <iostream>
#include "solver_pnp3d.cpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

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
	mwSize *dim0, *dim1, *dim2, *dim3;

	double *CNew;
	double *EzNew, *ErNew, *EqNew;
	double *JzNew, *JrNew, *JqNew;

	/* parse input arguments */

	C  = mxGetPr(prhs[0]);
	dim0 = (mwSize*)mxGetDimensions(prhs[0]);

	Ez = mxGetPr(prhs[1]);
	Er = mxGetPr(prhs[2]);
	Eq = mxGetPr(prhs[3]);

	Jz = mxGetPr(prhs[4]);
	dim1 = (mwSize*)mxGetDimensions(prhs[4]);
	Jr = mxGetPr(prhs[5]);
	dim2 = (mwSize*)mxGetDimensions(prhs[5]);
	Jq = mxGetPr(prhs[6]);
	dim3 = (mwSize*)mxGetDimensions(prhs[6]);


	m0 = (int)dim0[0];
 	n0 = (int)dim0[1];
 	l0 = (int)dim0[2];
 	m1 = (int)dim1[0];
	n1 = (int)dim1[1];
	l1 = (int)dim0[2];
 	m2 = (int)dim2[0];
	n2 = (int)dim2[1];
	l2 = (int)dim2[2];
	m3 = (int)dim3[0];
	n3 = (int)dim3[1];
	l3 = (int)dim3[2];

	// cout << dim0[0] << " " << dim0[1] << " " << dim0[2] << endl;
	// cout << dim1[0] << " " << dim1[1] << " " << dim1[2] << endl;
	// cout << dim2[0] << " " << dim2[1] << " " << dim2[2] << endl;

	z   = mxGetPr(prhs[7]);
	d_m = mxGetPr(prhs[8]);


  dqq = mxGetPr(prhs[9]);
	dx = mxGetPr(prhs[10]);
	dt = mxGetPr(prhs[11]);

	R = mxGetPr(prhs[12]);

	a = mxGetPr(prhs[13]);
	b = mxGetPr(prhs[14]);
	c = mxGetPr(prhs[15]);

	bulk = mxGetPr(prhs[16]);

	amo_z = mxGetPr(prhs[17]);
	amo_r = mxGetPr(prhs[18]);

	slp = mxGetPr(prhs[19]);

	n_display  = mxGetPr(prhs[20]);



	/* create output arguments */
    plhs[0] = mxCreateNumericArray(4, dim0, mxDOUBLE_CLASS, mxREAL);
    CNew = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(3, dim1, mxDOUBLE_CLASS, mxREAL);
    EzNew = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(3, dim2, mxDOUBLE_CLASS, mxREAL);
    ErNew = mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericArray(3, dim3, mxDOUBLE_CLASS, mxREAL);
    EqNew = mxGetPr(plhs[3]);
    plhs[4] = mxCreateNumericArray(4, dim1, mxDOUBLE_CLASS, mxREAL);
    JzNew = mxGetPr(plhs[4]);
    plhs[5] = mxCreateNumericArray(4, dim2, mxDOUBLE_CLASS, mxREAL);
    JrNew = mxGetPr(plhs[5]);
    plhs[6] = mxCreateNumericArray(4, dim3, mxDOUBLE_CLASS, mxREAL);
    JqNew = mxGetPr(plhs[6]);

		solver_pnp3d(C, Ez, Er, Eq, Jz, Jr, Jq,
			m0, n0, l0, m1, n1, l1, m2, n2, l2, m3, n3, l3, z, d_m,
			dqq, dx, dt, R, a, b, c, bulk, amo_z, amo_r, slp, n_display );

		/* copy results */

		for (int u = 0; u < (m0+1); u++ ) {
			for (int v = 0; v < n0; v++ ) {
				for (int w = 0; w < (l0+1); w++) {
					for (int k = 0; k < 3; k++) {
						// copy C & Jr
						if (u < m0 && w < l0) {
							CNew[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ]
								= C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ];
							JrNew[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ]
								= Jr[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ];
						}
						// copy Jz
						if (w < l0) {
							JzNew[ u + v*m1 + w*m1*n1 + k*m1*n1*l1 ]
								= Jz[ u + v*m1 + w*m1*n1 + k*m1*n1*l1 ];
						}
						// copy Jq
						if (u < m0) {
							JqNew[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ]
								= Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ];
						}
					} // end for-k
					// copy Ez
					if (w < l0) {
						EzNew[ u + v*m1 + w*m1*n1 ]
							= Ez[ u + v*m1 + w*m1*n1 ];
					}
					// copy Er
					if (u < m0 && w < l0) {
						ErNew[ u + v*m2 + w*m2*n2 ]
							= Er[ u + v*m2 + w*m2*n2 ];
					}
					// copy Eq
					if (u < m0) {
						EqNew[ u + v*m3 + w*m3*n3 ]
							= Eq[ u + v*m3 + w*m3*n3 ];
					}
				} // end for-w
			} // end for-v
		} // end for-u

}
