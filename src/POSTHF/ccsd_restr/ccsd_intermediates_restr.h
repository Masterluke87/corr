#ifndef CCSD_INTERMEDIATES_RESTR_H
#define CCSD_INTERMEDIATES_RESTR_H
#include "../ccsd.h"


struct cc_intermediates_restr
{
    double* Hik;
    double* Hca;
    double* Hck;

    double* Gik;
    double* Gca;

    double* Aijkl;
    double* Bcdab;
    double* Jicak;
    double* Kicka;

    double* tau;
    double* Pijab;
 };

void build_tau(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);
void build_Hik(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);
void build_Gik(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);

void build_Hca(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);
void build_Gca(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);

void build_Hck(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);

void build_Aikl(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);
void build_Bcdab(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);
void build_Jicak(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);
void build_Kicka(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);

void build_Pijab(cc_helper* CC,cc_intermediates_restr* CC_int,pHF* postHF);





#endif
