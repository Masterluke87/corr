double calc_mo2int_bf(int i, int j, int k, int l, int nroao, double* MOs,
		      long long int nrofint, long long int* sortcount, double* intval,
		      unsigned short* intnums);
double calc_1p_op_cis_od(int a, int b, double* MOs, int nroao, double* mat,
			 double* tmpvec);
double calc_1p_op_cis_d(int i, int f, double* MOs, int nroao, double* opmat,
			double* Pmat, int nroe);
void calc_mu_mat_cis(int cis_size, int nroao, int omo, int umo, int llim, int nroe,
		     double *cismat, double* mumat, double mucore,
		     double* MOs, double* Pmat, double* tmpvec_ao);
double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
		 long long int nrofint, long long int* sortcount,  double* intval,
		 unsigned short* intnums, std::ostream* outf);

//preints
int calc_mo2el_ind_x(int i, int j, int omo, int umo, int llim);
int calc_mo2el_ind_y(int k, int l, int omo, int umo, int llim);
double* precalc_ints_sd(int omo, int umo, int llim,
			int nroao, double* MOs,
			long long int nrofint, long long int* sortcount,  double* intval,
			unsigned short* intnums, std::ofstream* outf);
double get_precalc_ints_sd(int i, int j, int k, int l,
			   int omo, int umo, int llim, double* prec_ints);
double get_precalc_ints_fast(int i, int j, int k, int l,
			     int omo, int umo, int llim, double* prec_ints);
