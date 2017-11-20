

void run_scf(int nroao,int nroe, double* C,double* Pmat,double* Hmat,double* Fmat,
            unsigned short* intnums,double* intval, long long int* sortcount,
            long long int nrofint, double* Som12,int maxiter,double ion_rep);

void form_core_guess(int nroao,double* Fmat,double* Hmat,double * Som12,double* MOs,double* MOens);
