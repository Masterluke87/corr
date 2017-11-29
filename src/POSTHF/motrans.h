double calc_mo2int_bf(int i, int j, int k, int l, int nroao, double* MOs,
		      long long int nrofint, long long int* sortcount, double* intval,
		      unsigned short* intnums);
double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
		 long long int nrofint, long long int* sortcount,  double* intval,
		 unsigned short* intnums, std::ostream* outf);
