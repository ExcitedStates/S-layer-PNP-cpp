#include <mex.h>
#include <matrix.h>
#include <iostream>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	double *C;
	double *Ex, *Ey;
	double *Jx, *Jy;
	
	double *z, *d_m;

	double *dx, *dt;
	double *a, *b, *c;

	double *bulk;
	double *amo_y, *amo_x;
	double *slp;

	double *x_sym;
	double *n_display;

	double *no_slayer;
	
	int m0, n0, m1, n1, m2, n2;
	mwSize *dim0, *dim1, *dim2;

	double *CNew;
	double *ExNew, *EyNew;
	double *JxNew, *JyNew;

	/* parse input arguments */

	C  = mxGetPr(prhs[0]);
	dim0 = (mwSize*)mxGetDimensions(prhs[0]);

	Ex = mxGetPr(prhs[1]);
	Ey  = mxGetPr(prhs[2]);
	
	Jx = mxGetPr(prhs[3]);
	dim1 = (mwSize*)mxGetDimensions(prhs[3]);
	Jy = mxGetPr(prhs[4]);
	dim2 = (mwSize*)mxGetDimensions(prhs[4]);

	m0 = (int)dim0[0];
   	n0 = (int)dim0[1];
   	m1 = (int)dim1[0];
	n1 = (int)dim1[1];
   	m2 = (int)dim2[0];
	n2 = (int)dim2[1];

	// cout << dim0[0] << " " << dim0[1] << " " << dim0[2] << endl;
	// cout << dim1[0] << " " << dim1[1] << " " << dim1[2] << endl;
	// cout << dim2[0] << " " << dim2[1] << " " << dim2[2] << endl;


	z   = mxGetPr(prhs[5]);
	d_m = mxGetPr(prhs[6]);

	dx = mxGetPr(prhs[7]);
	dt = mxGetPr(prhs[8]);

	a = mxGetPr(prhs[ 9]);
	b = mxGetPr(prhs[10]);
	c = mxGetPr(prhs[11]);

	bulk = mxGetPr(prhs[12]);

	amo_y = mxGetPr(prhs[13]);
	amo_x = mxGetPr(prhs[14]);

	slp = mxGetPr(prhs[15]);

	x_sym      = mxGetPr(prhs[16]);
	n_display  = mxGetPr(prhs[17]);

	no_slayer  = mxGetPr(prhs[18]);

	/* create output arguments */
    plhs[0] = mxCreateNumericArray(3, dim0, mxDOUBLE_CLASS, mxREAL);
    CNew = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(2, dim1, mxDOUBLE_CLASS, mxREAL);
    ExNew = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(2, dim2, mxDOUBLE_CLASS, mxREAL);
    EyNew = mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericArray(3, dim1, mxDOUBLE_CLASS, mxREAL);
    JxNew = mxGetPr(plhs[3]);
    plhs[4] = mxCreateNumericArray(3, dim2, mxDOUBLE_CLASS, mxREAL);
    JyNew = mxGetPr(plhs[4]);




    for (int t = 0; t < (int)*n_display; t++) {


    	/* update Ex */

    	for (int u = 0; u < m1; u++) {
    		for (int v = 0; v < (int)*x_sym; v++) {
    			for (int w = 0; w < 3; w++) {
    				Ex[ u + v*m1 ] -=  *a * z[w] * Jx[ u + v*m1 + w*m1*n1 ]; 
    			}
    		}
    	}

    	/* update Ey */

    	for (int u = 0; u < m2; u++) {
    		for (int v = 0; v < (int)*x_sym; v++) {
    			for (int w = 0; w < 3; w++) {
    				Ey[ u + v*m2 ] -=  *a * z[w] * Jy[ u + v*m2 + w*m2*n2 ]; 
    			}
    		}
    	}

    	/* update Jx */

    	for (int u = 0; u < m1; u++) {
    		for (int v = 1; v < (int)*x_sym; v++) {
    			for (int w = 0; w < 3; w++) {
    				Jx[ u + v*m1 + w*m1*n1] = 0.5 * *b * d_m[w] * z[w] 
    				  * ( C[ u+v*m0+w*m0*n0 ] + C[ u+(v-1)*m0+w*m0*n0 ] ) * Ex[ u + v*m1 ];
    				Jx[ u + v*m1 + w*m1*n1] -= d_m[w] / *dx * ( C[ u+v*m0+w*m0*n0 ] - C[ u+(v-1)*m0+w*m0*n0 ] );
    			}
    		}
    	}

    	/* update Jy */

    	for (int u = 1; u < (m2-1); u++) {
    		for (int v = 1; v < (int)*x_sym; v++) {
    			for (int w = 0; w < 3; w++) {
    				Jy[ u + v*m2 + w*m2*n2] = 0.5 * *b * d_m[w] * z[w] 
    				  * ( C[ u+v*m0+w*m0*n0 ] + C[ (u-1)+v*m0+w*m0*n0 ] ) * Ey[ u + v*m2 ];
    				Jy[ u + v*m2 + w*m2*n2] -= d_m[w] / *dx * ( C[ u+v*m0+w*m0*n0 ] - C[ (u-1)+v*m0+w*m0*n0 ] );
    			}
    		}
    	}


    	/* diffusive BCs */

    	for (int v = 0; v < (int)*x_sym; v++) {
    		for (int w = 0; w < 3; w++) {
    			Jy[ 0 + v*m2 + w*m2*n2 ] = -d_m[w] / *dx * ( C[ 0 + v*m0 + w*m0*n0 ] - bulk[ w + v*3 ] ); 
    		}
    	}

    	/* Neumann BCs */

    	// cout << slp[0] << " " << slp[1] << " " << slp[2] << " " << slp[3] << endl;

	if (*no_slayer < 1.0) {
    	for (int w = 0; w < 3; w++) {

    		// vertical walls
    		for (int u = ((int)slp[0]-1); u < (int)slp[1]; u++) {
    			Jx[ u + (int)slp[2]*m1 + w*m1*n1 ] = 0;
    			Jx[ u + ((int)slp[3]-1)*m1 + w*m1*n1 ] = 0;
    		}

    		// horizontal surfaces
    		for (int v = 0; v < (int)slp[2]; v++) {
    			Jy[ ((int)slp[0]-1) + v*m2 + w*m2*n2 ] = 0;
    			Jy[ (int)slp[1] + v*m2 + w*m2*n2 ] = 0;
    		}
    		for (int v = ((int)slp[3]-1); v < n2; v++) {
    			Jy[ ((int)slp[0]-1) + v*m2 + w*m2*n2 ] = 0;
    			Jy[ (int)slp[1] + v*m2 + w*m2*n2 ] = 0;
    		}

    	}
	}

    	/* update C */

    	// Fake Jx at BC 
	    for (int u = 0; u < m0; u++) {
	    	for (int w = 0; w < 3; w++) {
	    		Jx[ u + (int)*x_sym * m1 + w*m1*n1 ]
	    		  = -Jx[ u + ((int)*x_sym - 1) * m1 + w*m1*n1 ];
	    	}
	    }

    	for (int u = 0; u < m0; u++) {
    		for (int v = 0; v < (int)*x_sym; v++) {
    			for (int w = 0; w < 3; w++) {
    				C[ u + v*m0 + w*m0*n0 ]
    				  -= *dt / *dx 
    				  * ( Jx[ u + (v+1)*m1 + w*m1*n1 ] - Jx[ u + v*m1 + w*m1*n1 ]
    				  + Jy[ (u+1) + v*m2 + w*m2*n2 ] - Jy[ u + v*m2 + w*m2*n2 ] );
    			}
    		}
    	}

    	/* reaction */

    	C[ ((int)*amo_y-1) + ((int)*amo_x-1) * m0 ] -= *c * C[ ((int)*amo_y-1) + ((int)*amo_x-1) * m0 ];

    }

	
    /* copy results */

    for (int u = 0; u < (m0+1); u++ ) {
    	for (int v = 0; v < (int)*x_sym; v++ ) {
    		for (int w = 0; w < 3; w++) {
    			// copy C & Jx
    			if (u < m0) {
    				CNew[ u + v*m0 + w*m0*n0 ] = C[ u + v*m0 + w*m0*n0 ];
    				JxNew[ u + v*m1 + w*m1*n1 ] = Jx[ u + v*m1 + w*m1*n1 ];
    			}
    			// copy Jy
    			JyNew[ u + v*m2 + w*m2*n2 ] = Jy[ u + v*m2 + w*m2*n2 ];
    		}
    		// copy Ex
    		if (u < m0) {
				ExNew[ u + v*m1 ] = Ex[ u + v*m1 ];
			}
			// copy Ey
			EyNew[ u + v*m2 ] = Ey[ u + v*m2 ];
    	}
    }

    /* free memory */
		
}
