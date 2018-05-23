void readBinHeader( string filename,
	int* t_offset,
	int* m0, int* n0, int* l0,
	int* m1, int* n1, int* l1,
	int* m2, int* n2, int* l2,
	int* m3, int* n3, int* l3
	) {

		int mode = -1;
		FILE *file;
		char filename_char[100];
		strcpy(filename_char, filename.c_str());
		file = fopen(filename_char , "rb");

		/* read header */
		fread(&mode, sizeof(int), 1, file);
		fread(t_offset, sizeof(int), 1, file);
		cout << "mode: " << mode << endl;
		cout << "time step: " << *t_offset << endl;

		/* read matrix shapes */
		fread(m0, sizeof(int), 1, file);
		fread(n0, sizeof(int), 1, file);
		fread(l0, sizeof(int), 1, file);
		fread(m1, sizeof(int), 1, file);
		fread(n1, sizeof(int), 1, file);
		fread(l1, sizeof(int), 1, file);
		fread(m2, sizeof(int), 1, file);
		fread(n2, sizeof(int), 1, file);
		fread(l2, sizeof(int), 1, file);
		fread(m3, sizeof(int), 1, file);
		fread(n3, sizeof(int), 1, file);
		fread(l3, sizeof(int), 1, file);

		fclose(file);

}

void readBin( string filename, const int n_i,
	double *C, double *Ez, double *Er, double *Eq,
	double *Jz, double *Jr, double *Jq,
	int* m0, int* n0, int* l0,
	int* m1, int* n1, int* l1,
	int* m2, int* n2, int* l2,
	int* m3, int* n3, int* l3,
	double *z, double *d_m,
	double *dqq, double *dx, double *dt, double *R,
	double *a, double *b, double *c,
	double *bulk, double *amo_z, double *amo_r, double *slp, double *n_display ) {

		int mode = -1;
		int t_offset = -1;
		FILE *file;
		char filename_char[100];
		strcpy(filename_char, filename.c_str());
		file = fopen(filename_char , "rb");
		cout << "opened " << filename_char << endl;

		/* read header */
		fread(&mode, sizeof(int), 1, file);
		fread(&t_offset, sizeof(int), 1, file);

		/* read matrix shapes */
		fread(m0, sizeof(int), 1, file);
		fread(n0, sizeof(int), 1, file);
		fread(l0, sizeof(int), 1, file);
		fread(m1, sizeof(int), 1, file);
		fread(n1, sizeof(int), 1, file);
		fread(l1, sizeof(int), 1, file);
		fread(m2, sizeof(int), 1, file);
		fread(n2, sizeof(int), 1, file);
		fread(l2, sizeof(int), 1, file);
		fread(m3, sizeof(int), 1, file);
		fread(n3, sizeof(int), 1, file);
		fread(l3, sizeof(int), 1, file);

		/* read C profile */
		int size0 = 0;
		size0 = *m0 * *n0 * *l0 * n_i;
		// cout << "C has " << size0 << " elements." << endl;
		fread(C, sizeof(double), size0, file);

		/* read E profiles */
		int size1 = 0;
		size1 = *m1 * *n1 * *l1;
		cout << "Ez has " << size1 << " elements." << endl;
		fread(Ez, sizeof(double), size1, file);
		int size2 = 0;
		size2 = *m2 * *n2 * *l2;
		cout << "Er has " << size2 << " elements." << endl;
		fread(Er, sizeof(double), size2, file);
		int size3 = 0;
		size3 = *m3 * *n3 * *l3;
		cout << "Eq has " << size3 << " elements." << endl;
		fread(Eq, sizeof(double), size3, file);

		/* read J profiles */
		size1 = *m1 * *n1 * *l1 * n_i;
		cout << "Jz has " << size1 << " elements." << endl;
		fread(Jz, sizeof(double), size1, file);
		size2 = *m2 * *n2 * *l2 * n_i;
		cout << "Jr has " << size2 << " elements." << endl;
		fread(Jr, sizeof(double), size2, file);
		size3 = *m3 * *n3 * *l3 * n_i;
		cout << "Jq has " << size3 << " elements." << endl;
		fread(Jq, sizeof(double), size3, file);

		/* read constants */
		fread(z, sizeof(double), n_i, file);
		fread(d_m, sizeof(double), n_i, file);
		fread(dqq, sizeof(double), 1, file);
		fread(dx, sizeof(double), 1, file);
		fread(dt, sizeof(double), 1, file);

		/* read R */
		size0 = *m0 * *n0 * (*l0+1);
		cout << "R has " << size0 << " elements." << endl;
		fread(R, sizeof(double), size0, file);

		/* read coefficients */
		fread(a, sizeof(double), 1, file);
		fread(b, sizeof(double), 1, file);
		fread(c, sizeof(double), 1, file);

		/* read bulk */
		size0 =  *n0 * *l0 * n_i;
		cout << "bulk has " << size0 << " elements." << endl;
		fread(bulk, sizeof(double), size0, file);

		/* read misc. */
		fread(amo_z, sizeof(double), 1, file);
		fread(amo_r, sizeof(double), 1, file);
		fread(slp, sizeof(double), 3, file);
		fread(n_display, sizeof(double), 1, file);

		fclose(file);

}

void writeBin( string filename, const int n_i, int t_offset,
	double *C, double *Ez, double *Er, double *Eq,
	double *Jz, double *Jr, double *Jq,
	int* m0, int* n0, int* l0,
	int* m1, int* n1, int* l1,
	int* m2, int* n2, int* l2,
	int* m3, int* n3, int* l3
	) {

		int mode = 0;
		int t = -1;
		FILE *file;
		char filename_char[100];
		strcpy(filename_char, filename.c_str());
		file = fopen(filename_char , "wb");

		/* write header */
		fwrite(&mode, sizeof(int), 1, file);
		fwrite(&t_offset, sizeof(int), 1, file);

		/* write matrix shapes */
		fwrite(m0, sizeof(int), 1, file);
		fwrite(n0, sizeof(int), 1, file);
		fwrite(l0, sizeof(int), 1, file);
		fwrite(m1, sizeof(int), 1, file);
		fwrite(n1, sizeof(int), 1, file);
		fwrite(l1, sizeof(int), 1, file);
		fwrite(m2, sizeof(int), 1, file);
		fwrite(n2, sizeof(int), 1, file);
		fwrite(l2, sizeof(int), 1, file);
		fwrite(m3, sizeof(int), 1, file);
		fwrite(n3, sizeof(int), 1, file);
		fwrite(l3, sizeof(int), 1, file);

		/* write C profile */
		int size0 = 0;
		size0 = *m0 * *n0 * *l0 * n_i;
		// cout << "C has " << size0 << " elements." << endl;
		fwrite(C, sizeof(double), size0, file);

		/* write E profiles */
		int size1 = 0;
		size1 = *m1 * *n1 * *l1;
		// cout << "Ez has " << size1 << " elements." << endl;
		fwrite(Ez, sizeof(double), size1, file);
		int size2 = 0;
		size2 = *m2 * *n2 * *l2;
		// cout << "Er has " << size2 << " elements." << endl;
		fwrite(Er, sizeof(double), size2, file);
		int size3 = 0;
		size3 = *m3 * *n3 * *l3;
		// cout << "Eq has " << size3 << " elements." << endl;
		fwrite(Eq, sizeof(double), size3, file);

		/* write J profiles */
		size1 = *m1 * *n1 * *l1 * n_i;
		// cout << "Jz has " << size1 << " elements." << endl;
		fwrite(Jz, sizeof(double), size1, file);
		size2 = *m2 * *n2 * *l2 * n_i;
		// cout << "Jr has " << size2 << " elements." << endl;
		fwrite(Jr, sizeof(double), size2, file);
		size3 = *m3 * *n3 * *l3 * n_i;
		// cout << "Jq has " << size3 << " elements." << endl;
		fwrite(Jq, sizeof(double), size3, file);

		fclose(file);

}
