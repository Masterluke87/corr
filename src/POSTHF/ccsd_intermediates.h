#ifndef CCSD_INTERMEDIATES_H
#define CCSD_INTERMEDIATES_H

#include "ccsd.h"

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

    double* pMem;
};


struct cc_intermediates_restr
{
    double* Hik;
    double* Hca;
    double* Hck;
};

void allocate_intermediates_memory(cc_helper *CC,cc_intermediates* CC_int);
void build_Fae(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF);
void build_Fmi(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF);
void build_Fme(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF);
void build_Wmnij(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF);
void build_Wabef(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF);
void build_Wmbej(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF);





#endif
