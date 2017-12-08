#include <iostream>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include "ops_rhf.h"
#include <cblas.h>
#define PWIDTH_L 12
#define PWIDTH_R 16

extern "C"
{
int dgemm_(char *, char *,
           int *, int *, int *,
           double *, double *,
           int *,double *,
           int *,double *, double *, int *);
int daxpy_(int*,double*,double*,int*,double*, int*);
}

void form_core_guess(int nroao,double* Fmat, double* Hmat,double * Som12,double* MOs,double* MOens)
{
    for (int i =0; i<nroao*nroao; i++)
        Fmat[i] = Hmat[i];
    double *tmpmat = new double[nroao*nroao];

    diag_Fmat(nroao, Fmat,MOs,MOens,Som12, tmpmat);
    delete[] tmpmat;
}








void build_Pmat(double *C,int nroao, int nroe, double *Pmat,double damp)
{
    double alpha=2.0 * (1-damp);
    int nocc = nroe/2;
    double c = 0.0;
    double PmatOld[nroao*nroao];

    for (int i = 0; i < nroao*nroao; i++) {
        PmatOld[i] = Pmat[i];
    }
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

double twoE(int nroao,double* Pmat,double* Fmat){
    double twoE = 0.0;
    for (size_t mu = 0; mu < nroao; mu++) {
        for (size_t nu = 0; nu < nroao; nu++) {
            twoE += Pmat[mu*nroao + nu] * Fmat[mu*nroao + nu];
        }
    }
    return twoE;
}



double calc_diis_error(int nroao,double* Fmat,double* Pmat,double* Smat){

    double* SDF = new double[nroao*nroao];
    double* FDS = new double[nroao*nroao];
    memset(SDF,0,nroao*nroao*sizeof(double));
    memset(FDS,0,nroao*nroao*sizeof(double));
    for(int i=0; i<nroao; i++) {
        for(int j=0; j<nroao; j++) {
            for (int mu = 0; mu < nroao; mu++) {
                for (int nu = 0; nu < nroao; nu++) {
                    SDF[i*nroao+j] += Smat[i*nroao+nu] *Pmat[nu*nroao+mu] *  Fmat[mu*nroao+j];
                    FDS[i*nroao+j] += Fmat[i*nroao+nu] *Pmat[nu*nroao+mu] *  Smat[mu*nroao+j];

                }
            }

        }
    }
    double error = 0.0;
    for (int i=0;i<nroao*nroao;i++){
        error += std::fabs(SDF[i] - FDS[i]);
    }
    delete[] SDF;
    delete[] FDS;
    return error/(nroao*nroao);
}
void run_scf(int nroao,int nroe, double* C, double* Pmat,double* Hmat,double* Smat,double* Fmat,double* MOens,
             unsigned short* intnums, double* intvals, long long int* sortcount,
             long long int nrofint, double* Som12, int maxiter,double ion_rep){
    std::cout << "+++++++++++SCF++++++++++" << '\n';


    double* tmpmat = new double[nroao*nroao];

    double Escf  = 0.0;
    double DE    = 10.0;
    double Eold  = 100.0;
    double damp  = 0.5;
    double start = 0.0;
    double end   = 0.0;
    double error = 0.0;
    int iter=0;
    build_Pmat(C,nroao,nroe,Pmat,0.0);
    std::cout<<std::fixed<<std::setprecision(10);
    std::cout<<std::setw(-4)<<"Iter"<<std::setw( 16 )<<"ESCF" << std::setw( 16 )<<"DE"<<std::setw( 16 )<<"ERROR"<<std::setw( 16 )<<"t[s]"<<"\n";
    while (iter  < maxiter && std::fabs(DE) > 1E-8) {
        start = omp_get_wtime();
        //with guess
        build_Pmat(C,nroao,nroe,Pmat,damp);
        build_Fmat(nroao,Fmat,Pmat,Hmat,intvals,intnums,sortcount,nrofint);
        error = calc_diis_error(nroao,Fmat,Pmat,Smat);
        Escf = Calc_e_el(nroao,Fmat,Pmat,Hmat);
        DE = Escf - Eold;
        Eold = Escf;
        end = omp_get_wtime();
        std::cout<<std::setw(4)<<iter<<":"<<
            std::setw( 16 )<< Escf+ion_rep<<std::setw( 16 )<<DE <<std::setw( 16 )<<error<<std::setw( 16 )<<(end-start)<<'\n';
        diag_Fmat(nroao, Fmat,C,MOens,Som12, tmpmat);


        // build_Pmat(C,nroao,nroe,Pmat,damp);
        // build_Fmat(nroao,Fmat,Pmat,Hmat,intvals,intnums,sortcount,nrofint);
        // diag_Fmat(nroao, Fmat,C,MOens,Som12, tmpmat);
        iter++;
    }
    build_Pmat(C,nroao,nroe,Pmat,damp);
    build_Fmat(nroao,Fmat,Pmat,Hmat,intvals,intnums,sortcount,nrofint);


    std::cout<<"\nHF-Results:\n";
    std::cout<<"-----------\n";
    double one,two;
    one = oneE(nroao,Pmat,Hmat);
    two = twoE(nroao,Pmat,Fmat);
    std::cout <<std::setw( PWIDTH_L ) << "  ONE E:" <<std::setw( PWIDTH_R ) <<one+ion_rep << '\n';
    std::cout <<std::setw( PWIDTH_L ) << "  TWO E:" <<std::setw( PWIDTH_R )<<(two-one)/2 << '\n';
    std::cout <<std::setw( PWIDTH_L ) << "TOTAL E:" <<std::setw( PWIDTH_R )<<(one+two)/2+ion_rep<< '\n';

    delete[] tmpmat;
    //delete MOens;
}
