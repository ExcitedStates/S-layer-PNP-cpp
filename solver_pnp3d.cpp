void solver_pnp3d(double *C, double *Ez, double *Er, double *Eq,
	double *Jz, double *Jr, double *Jq,
	int m0, int n0, int l0, int m1, int n1, int l1,
	int m2, int n2, int l2, int m3, int n3, int l3,
	double *z, double *d_m, double *dqq, double *dx, double *dt, double *R,
	double *a, double *b, double *c, double *bulk,
	double *amo_z, double *amo_r, double *slp, double *n_display
	) {

	for (int t = 0; t < (int)*n_display; t++) {

		/* update Ez */

		for (int u = 0; u < m1; u++) {
			for (int v = 0; v < n1; v++) {
				for (int w = 0; w < l1; w++) {
					for (int k = 0; k < 3; k++) {
						Ez[ u + v*m1 + w*m1*n1 ] -=  *a * z[k] * Jz[ u + v*m1 + w*m1*n1 + k*m1*n1*l1 ];
					}
				}
			}
		}

		/* update Er */

		for (int u = 0; u < m2; u++) {
			for (int v = 0; v < n2; v++) {
				for (int w = 0; w < l2; w++) {
					for (int k = 0; k < 3; k++) {
						Er[ u + v*m2 + w*m2*n2 ] -=  *a * z[k] * Jr[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ];
					}
				}
			}
		}

		/* update Eq */

		for (int u = 0; u < m3; u++) {
			for (int v = 0; v < n3; v++) {
				for (int w = 0; w < l3; w++) {
					for (int k = 0; k < 3; k++) {
						Eq[ u + v*m2 + w*m2*n2 ] -=  *a * z[k] * Jq[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ];
					}
				}
			}
		}

		/* update Jz */

		for (int u = 1; u < (m1-1); u++) {
			for (int v = 0; v < n1; v++) {
				for (int w = 0; w < l1; w++) {
					for (int k = 0; k < 3; k++) {
						Jz[ u + v*m1 + w*m1*n1 + k*m1*n1*l1 ]
							= 0.5 * *b * d_m[k] * z[k]
	  				  * ( C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] + C[ (u-1) + v*m0 + w*m0*n0 + k*m0*n0*l0 ] )
							* Ez[ u + v*m1 + w*m1*n1 ];
	  				Jz[ u + v*m1 + w*m1*n1 + k*m1*n1*l1 ]
							-= d_m[k] / *dx
							* ( C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] - C[ (u-1) + v*m0 + w*m0*n0 + k*m0*n0*l0 ] );
					}
				}
			}
		}

		/* update Jr */

		for (int u = 0; u < m2; u++) {
			for (int v = 0; v < (n2-1); v++) {
				for (int w = 0; w < l2; w++) {
					for (int k = 0; k < 3; k++) {
						Jr[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ]
							= 0.5 * *b * d_m[k] * z[k]
	  				  * ( C[ u + (v+1)*m0 + w*m0*n0 + k*m0*n0*l0 ] + C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] )
							* Er[ u + v*m2 + w*m2*n2 ];
	  				Jr[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ]
							-= d_m[k] / *dx
							* ( C[ u + (v+1)*m0 + w*m0*n0 + k*m0*n0*l0 ] - C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] );
					}
				}
			}
		}

		/* update Jq */

		for (int u = 0; u < m3; u++) {
			for (int v = 0; v < n3; v++) {
				for (int k = 0; k < 3; k++) {
					for (int w = 1; w < (l3-1); w++) {
						Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ]
							= 0.5 * *b * d_m[k] * z[k]
	  				  * ( C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] + C[ u + v*m0 + (w-1)*m0*n0 + k*m0*n0*l0 ] )
							* Eq[ u + v*m3 + w*m3*n3 ];
	  				Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ]
							-= d_m[k] / *dx
							* ( C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] - C[ u + v*m0 + (w-1)*m0*n0 + k*m0*n0*l0 ] );
					} // end for-w
					// Jq BC-1
					int w = 0;
					Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ]
						= 0.5 * *b * d_m[k] * z[k]
						* ( C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] + C[ u + v*m0 + (l0-1)*m0*n0 + k*m0*n0*l0 ] )
						* Eq[ u + v*m3 + w*m3*n3 ];
					Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ]
						-= d_m[k] / *dx
						* ( C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ] - C[ u + v*m0 + (l0-1)*m0*n0 + k*m0*n0*l0 ] );
				  // Jq BC-2
					w = l3-1;
					Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ]
						= 0.5 * *b * d_m[k] * z[k]
						* ( C[ u + v*m0 + 0*m0*n0 + k*m0*n0*l0 ] + C[ u + v*m0 + (w-1)*m0*n0 + k*m0*n0*l0 ] )
						* Eq[ u + v*m3 + w*m3*n3 ];
					Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ]
						-= d_m[k] / *dx
						* ( C[ u + v*m0 + 0*m0*n0 + k*m0*n0*l0 ] - C[ u + v*m0 + (w-1)*m0*n0 + k*m0*n0*l0 ] );
				} // end for-k
			} // end for-v
		} // end for-u

		/* diffusive BCs */

		for (int v = 0; v < n1; v++) {
			for (int w = 0; w < l1; w++) {
				for (int k = 0; k < 3; k++) {
					Jz[ 0 + v*m1 + w*m1*n1 + k*m1*n1*l1 ]
						= -d_m[k] / *dx
						* ( C[ 0 + v*m0 + w*m0*n0 + k*m0*n0*l0 ] - bulk[ 0 + v*1 + w*1*n0 + k*1*n0*l0 ] );
				}
			}
		}

		/* Neumann BCs */

		// cout << slp[0] << " " << slp[1] << " " << slp[2] << " " << slp[3] << endl;

		for (int k = 0; k < 3; k++) {
	  	for (int w = 0; w < l1; w++) {

	  		// vertical walls
	  		for (int u = ((int)slp[0]-1); u < (int)slp[1]; u++) {
	  			Jr[ u + ((int)slp[2]-2)*m2 + w*m2*n2 + k*m2*n2*l2 ] = 0;
	  		}

	  		// horizontal surfaces
	  		for (int v = ((int)slp[2]-1); v < n1; v++) {
	  			Jz[ ((int)slp[0]-1) + v*m1 + w*m1*n1 + k*m1*n1*l1 ] = 0;
	  			Jz[ (int)slp[1] + v*m1 + w*m1*n1 + k*m1*n1*l1 ] = 0;
	  		}
	  	}
		}

		/* update C */

		for (int k = 0; k < 3; k++) {
			for (int u = 0; u < m0; u++) {
				for (int w = 0; w < l0; w++) {
					for (int v = 0; v < n0; v++) {
						// Jz
						C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ]
							-= *dt / *dx
							* ( Jz[ (u+1) + v*m1 + w*m1*n1 + k*m1*n1*l1 ] - Jz[ u + v*m1 + w*m1*n1 + k*m1*n1*l1 ] );
						if (v > 0) {
							// Jr-1
							C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ]
								-= *dt / *dx
								* ( Jr[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ] - Jr[ u + (v-1)*m2 + w*m2*n2 + k*m2*n2*l2 ] );
							// Jr-2
							C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ]
								-= 0.5 * *dt
								/ R[ u + v*m0 + (w+1)*m0*n0 ]
								* ( Jr[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ] + Jr[ u + (v-1)*m2 + w*m2*n2 + k*m2*n2*l2 ] );
							// Jq
							C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ]
								-= *dt / *dqq
								/ R[ u + v*m0 + (w+1)*m0*n0 ]
								* ( Jq[ u + v*m3 + (w+1)*m3*n3 + k*m3*n3*l3 ] - Jq[ u + v*m3 + w*m3*n3 + k*m3*n3*l3 ] );
						} else {
	                            C[ u + v*m0 + w*m0*n0 + k*m0*n0*l0 ]
	                                -= *dt / *dx
								* 2 * Jr[ u + v*m2 + w*m2*n2 + k*m2*n2*l2 ];
	                        }

					}
					// special case: r = 0
					C[ u + 0*m0 + w*m0*n0 + k*m0*n0*l0 ]
						-= 2 * *dt / *dx
						* Jr[ u + 0*m0 + w*m0*n0 + k*m0*n0*l0 ];
				}
			}
		}

		/* reaction */

		for (int w = 0; w < l0; w++) {
			C[ ((int)*amo_z-1) + ((int)*amo_r-1) * m0 + w*m0*n0 ]
				-= *c * C[ ((int)*amo_z-1) + ((int)*amo_r-1) * m0 + w*m0*n0 ];
	                            if (C[ ((int)*amo_z-1) + ((int)*amo_r-1) * m0 + w*m0*n0 ] <= 0) {
	                                    C[ ((int)*amo_z-1) + ((int)*amo_r-1) * m0 + w*m0*n0 ] = 0;
	                            }
		}

	} // end for-t

}
