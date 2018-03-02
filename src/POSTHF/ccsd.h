#include "../IO/ops_io.h"


void ccsd_ur(systeminfo* sysinfo,OEints* onemats,pHF* postHF);


struct cc_intermediates
{
    double* Fae;//    = new double [nvir*nvir];
    double* Fmi;//    = new double [nocc*nocc];
    double* Fme;//    = new double [nocc*nvir];

    double* Wmnij;//  = new double [nocc*nocc*nocc*nocc];
    double* Wabef;//  = new double [nvir*nvir*nvir*nvir];
    double* Wmbej;//  = new double [nocc*nvir*nvir*nocc];

    double* tau_s;//  = new double [nocc*nocc*nvir*nvir];
    double* tau;//    = new double [nocc*nocc*nvir*nvir];
};


struct cc_helper{

    int nocc;
    int norb;
    int nvir;
    int nroao;
    double ion_rep;

    double* Hmat;
    double* MOs; //
    double* f;   //  fock matrix so basis
    double* h;
    double* T2;  //   = new double[nocc*nocc*nvir*nvir];
    double* T2n; //   = new double[nocc*nocc*nvir*nvir];
    double* T1;  //   = new double[nocc*nvir];
    double* T1n; //   = new double[nocc*nvir];

};
