#include "ccsd.h"
#include "../UTIL/util.h"


void ccsd_ur_ri(systeminfo* sysinfo,OEints* onemats,pHF* postHF){
    cc_helper *CC = new cc_helper;
    CC->nocc  = sysinfo->nroe;
    CC->nroao = sysinfo->nroao;
    CC->nvir  = (2*(sysinfo->nroao) - sysinfo->nroe);
    CC->norb  = 2*(sysinfo->nroao);
    CC->MOs   = onemats->MOs;
    CC->Hmat  = onemats->Hmat;
    CC->ion_rep = sysinfo->ion_rep;
    CC->Ecc     = 0.0;
    CC->Ecc_old = 0.0;
    CC->DE      = 100000.0;
    CC->iter    = 0;

}
