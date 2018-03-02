#ifndef RI_H
#define RI_H


#include <libint2.hpp>
#include <iostream>
#include <cblas.h>
#include <lapacke.h>
#include "../IO/ops_io.h"

void calculate_ri(libint2::BasisSet &obs, libint2::BasisSet &dfbs, double *BPQ);
void transform_ri(systeminfo* sysinfo, OEints* onemats,pHF* postHF);

#endif
