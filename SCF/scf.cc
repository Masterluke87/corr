#include <iostream>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cblas.h>
#include <omp.h>

#include "ops_rhf.h"

extern "C"
 {
   int dgemm_(char *, char *,
            int *, int *, int *,
             double *, double *,
              int *,double *,
              int *,double *, double *, int *);
   int daxpy_(int* ,double*,double*,int*,double*, int*);
 }

void build_Pmat(double *C,int nroao, int nroe, double *Pmat,double damp)
{
  double alpha=2.0 * (1-damp);
  int nocc = nroe/2;
  double c = 0.0;
  double PmatOld[nroao*nroao];

  //for (int i = 0; i < nroao*nroao; i++) {
  //  PmatOld[i] = Pmat[i];
  //}
  std::memcpy(PmatOld,Pmat,nroao*nroao*sizeof(double));

  dgemm_("N","T",&nroao,&nroao,&nocc,
        &alpha,C,&nroao,C,&nroao,&c,Pmat,&nroao);
  int size = nroao*nroao;
  double x = damp;
  int inc =1;
  daxpy_(&size,&x,PmatOld,&inc,Pmat,&inc);

  /*
  int mu,nu;
  for(int x = 0; x < nroao*nroao; x++) Pmat[x] = 0.;
  for (int i = 0; i < nroe/2; i++) {
    for ( mu = 0; mu < nroao; mu++) {
      for ( nu = 0; nu < nroao; nu++) {
        Pmat[mu*nroao + nu] += 2*C[i*nroao + mu ]* C[i*nroao + nu ];
      }
    }
  }
  double diff = 0.0;
  */
  //Blas checker
  //for (size_t x = 0; x < nroao*nroao; x++) {
  //  diff += std::fabs(Pmat[x] - PMat2[x]);
  //}
  //std::cout << diff << '\n';

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
  double MOens [nroao];
  double Escf = 0.0;
  double DE   = 10.0;
  double Eold = 0.0;
  double damp = 0.5;
  double start = 0.0;
  double end  = 0.0;
  int iter=0;
  std::cout<<std::fixed<<std::setprecision(10);
  std::cout<<std::setw(-4)<<"Iter"<<std::setw( 16 )<<"ESCF" << std::setw( 16 )<<"DE"<<std::setw( 16 )<<"t[s]"<<"\n";
  while (iter  < maxiter && std::fabs(DE)>1E-8) {
    start = omp_get_wtime();
    build_Pmat(C,nroao,nroe,Pmat,damp);
    diag_Fmat(nroao, Fmat,C,MOens,Som12, tmpmat);
    build_Fmat(nroao,Fmat,Pmat,Hmat,intvals,intnums,sortcount,nrofint);
    Escf = Calc_e_el(nroao,Fmat,Pmat,Hmat);
    DE = Escf - Eold;
    Eold = Escf;
    end = omp_get_wtime();
    std::cout<<std::setw(4)<<iter<<":"<<
    std::setw( 16 )<< Escf+ion_rep<<std::setw( 16 )<<DE <<std::setw( 16 )<<(end-start)<<'\n';
    iter++;
}


  delete tmpmat;
  //delete MOens;
}
