#ifndef CCSD_H
#define CCSD_H

#include "../IO/ops_io.h"

double get_integral(double* ints,long long int &istep,long long int &jstep,long long int &kstep,int  &i,int  &j,int &k,int &l);
void ccsd_ur   (systeminfo* sysinfo,OEints* onemats,pHF* postHF);
void ccsd_restr(systeminfo* sysinfo,OEints* onemats,pHF* postHF);

struct cc_helper{

    int nocc;
    int norb;
    int nvir;
    int nroao;
    double ion_rep;

    double* Hmat;
    double* MOs; //
    double* f;   //  fock matrix so basis
    double* f_ri;
    double* h;

    double* pMem;    //
    double* T2;     //   = new double[nocc*nocc*nvir*nvir];
    double* T2n;    //   = new double[nocc*nocc*nvir*nvir];
    double* T1;     //   = new double[nocc*nvir];
    double* T1n;    //   = new double[nocc*nvir];

    double Ecc;
    double Ecc_old;
    double DE;

    int iter;

};


#endif
