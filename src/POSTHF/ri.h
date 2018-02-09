#ifndef RI_H
#define RI_H


#include <libint2.hpp>
#include <iostream>
#include <cblas.h>
#include <lapacke.h>

void calculate_ri(libint2::BasisSet &obs, libint2::BasisSet &dfbs, double *BPQ);

#endif
