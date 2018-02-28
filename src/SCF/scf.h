#ifndef SCF_H
#define SCF_H


void run_scf(systeminfo *sysinfo,OEints *onemats,TEints* twomats);

void form_core_guess(int nroao,double* Fmat,double* Hmat,double * Som12,double* MOs,double* MOens);

#endif
