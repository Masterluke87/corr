#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "ops_rhf.h"


void build_Pmat(double *C,int nroao, int nroe, double *Pmat)
{
  int mu,nu;
  for(int x = 0; x < nroao*nroao; x++) Pmat[x] = 0.;
  for (int i = 0; i < nroe/2; i++) {
    for ( mu = 0; mu < nroao; mu++) {
      for ( nu = 0; nu < nroao; nu++) {
        Pmat[mu*nroao + nu] += C[i*nroao + mu ]* C[i*nroao + nu ];
      }
    }
  }
}

double oneE(int nroao, double* Pmat, double* Hmat){
  double oneE = 0.0;
    for (size_t mu = 0; mu < nroao; mu++) {
      for (size_t nu = 0; nu < nroao; nu++) {
        oneE += Pmat[mu*nroao + nu] * Hmat[mu*nroao + nu];
      }
    }
    return oneE;
}



void run_scf(int nroao,int nroe, double* C, double* Pmat,double* Hmat,double* Fmat,
            unsigned short* intnums, double* intvals, long long int* sortcount,
            long long int nrofint, double* Som12, int maxiter,double ion_rep){
  std::cout << "**SCF::run_scf**" << '\n';

  double* tmpmat = new double[nroao*nroao];
  double* Pmat_old = new double[nroao*nroao];
  double* MOens  = new double[nroao];
  double Escf = 0.0;
  double DE   = 10.0;
  double Eold = 0.0;
  double damp = 0.5;
  int iter=0;
  build_Pmat_dscf(nroao,nroe, Pmat, Pmat_old, C,  damp);
  build_Fmat(nroao,Fmat,Pmat,Hmat,intvals,intnums,sortcount,nrofint);
  diag_Fmat(nroao, Fmat,C,MOens,Som12, tmpmat);
  std::cout<<std::fixed<<std::setprecision(10);
  std::cout<<std::setw(-4)<<"Iter"<<std::setw( 16 )<<"ESCF" << std::setw( 16 )<<"DE\n";
  while (iter  < maxiter && std::fabs(DE)>1E-8) {
    std::cout << "oneE:" << oneE(nroao,Pmat,Hmat) <<'\n';
    build_Pmat_dscf(nroao,nroe, Pmat, Pmat_old, C,  damp);
    diag_Fmat(nroao, Fmat,C,MOens,Som12, tmpmat);
    build_Fmat(nroao,Fmat,Pmat,Hmat,intvals,intnums,sortcount,nrofint);
    Escf = Calc_e_el(nroao,Fmat,Pmat,Hmat);

    DE = Escf - Eold;

    std::cout<<std::setw(4)<<iter<<":"<<
    std::setw( 16 )<< Escf+ion_rep<<std::setw( 16 )<<DE << '\n';
    std::cout << "oneE:" << oneE(nroao,Pmat,Hmat) <<'\n';
    iter++;
}


  delete tmpmat;
  delete MOens;
  delete Pmat_old;
}
