#include "../IO/ops_io.h"



void MOtrans(systeminfo *sysinfo,OEints* onemats,TEints *twomats,pHF* postHF);
double calc_mo2int_bf(int i, int j, int k, int l, int nroao, double* MOs,
                      long long int nrofint, long long int* sortcount, double* intval,
                      unsigned short* intnums);
double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
                 long long int nrofint, long long int* sortcount,  double* intval,
                 unsigned short* intnums, std::ostream* outf);
void read_transform_ri(std::string prefix,  //prefix to find the file
                       int nroe,            //nr of electrons
                       int nroao,           //nr of bsf in aobasis
                       int naux_2, double *MOs,          //nr of aux basis functions
                       double* Bia                 //output - transformed array containing mointegrals
                       );



void MOtrans(systeminfo *sysinfo,OEints* onemats,TEints *twomats,pHF* postHF);
void build_FMo(systeminfo *sysinfo,OEints* onemats, pHF* postHF);
