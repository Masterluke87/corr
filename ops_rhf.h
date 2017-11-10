double calc_r_ab(int a, int b, double* coord);
double calc_ion_rep(int nroa, double* coord, double* charges);
void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
void calc_mu_core(int nroa, double* coord, double* charges, double* point,
		  double* mu_core);
double calc_S12(int nroao, double* Smat, double* Som12, double* tmpmat,
		double* tmpvecs, double* tmpvals);
double calc_Sp12(int nroao, double* Smat, double* Sop12, double* tmpmat,
		 double* tmpvecs, double* tmpvals);
void   diag_Fmat(int nroao, double* Fmat, double* MOs,
		 double* MOens, double* Som12, double* tmpmat);
void   build_Fmat(int nroao, double* Fmat, double* Pmat, double* Hmat,
		  double* intval, unsigned short* intnums,
		  long long int* sortcount, long long int nrofint);
double Calc_e_el(int nroao, double* Fmat, double* Pmat, double* Hmat);
double calc_op_1el(int nroao, double* opmat, double* Pmat);
