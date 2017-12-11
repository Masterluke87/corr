#include <iostream>
#include <cstring>
#include <memory>
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

double dnrm2_(int *,double*, int* );

void dgels_( char* , int* , int* , int* , double* , int* ,
                double* , int* , double* , int* , int* );

}

extern "C" void  dgesv_(int* n, int*  nrhs,double* a,int* lda,int *ipiv,double* b,int* ldb,int* info);


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

    double* SDF  = new double[nroao*nroao];
    double* FDS  = new double[nroao*nroao];
    double* grad = new double[nroao*nroao];
    memset(SDF,0,nroao*nroao*sizeof(double));
    memset(FDS,0,nroao*nroao*sizeof(double));
    memset(grad,0,nroao*nroao*sizeof(double));

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
        grad[i] = (SDF[i] - FDS[i]);
    }
    delete[] SDF;
    delete[] FDS;
    int ndim = nroao*nroao;
    int incx = 1;
    double norm = 0.0;
    norm = dnrm2_(&ndim,grad,&incx);
    delete[] grad;
    return norm;
}
void run_scf(int nroao,int nroe, double* C, double* Pmat,double* Hmat,double* Smat,double* Fmat,double* MOens,
             unsigned short* intnums, double* intvals, long long int* sortcount,
             long long int nrofint, double* Som12, int maxiter,double ion_rep){
    std::cout << "+++++++++++SCF++++++++++" << '\n';


    double* tmpmat = new double[nroao*nroao];

    double diis_thresh = 0.01;
    bool   diis        = false;
    int    diis_size   = 4;
    double* diis_mem   = new double[2*diis_size*nroao*nroao];
    memset(diis_mem,0,diis_size*nroao*nroao*sizeof(double));



    double Escf  = 0.0;
    double DE    = 10.0;
    double Eold  = 100.0;
    double damp  = 0.5;
    double start = 0.0;
    double end   = 0.0;

    std::vector<double> diis_error(diis_size,0.0);
    std::vector<double*> diis_mat(diis_size);

    for (int i=0;i<diis_size;i++){
        diis_mat[i] = &diis_mem[i*nroao*nroao];
    }


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


        for(int i=1;i<diis_size;i++)
              diis_error[i-1] = diis_error[i];

        for(int i=1;i<diis_size;i++){
            for(int j=0;j<nroao*nroao;j++){
                diis_mat[i-1][j] = diis_mat[i][j];
            }
        }

        for(int j=0;j<nroao*nroao;j++)
            diis_mat.back()[j] = Fmat[j];



        diis_error.back() = error;


        if (error<diis_thresh){
            if (!diis){
                std::cout<<"Starting DIIS \n";
                diis = true;
            }
            damp =0.0;
            std::cout<<"E: [";
            for(auto a:diis_error)
                std::cout<<a<<", ";
            int d_size = (diis_size+1);
            double* B = new double[d_size*d_size];
            double* c = new double[d_size];

            for (int i=0;i<diis_size;i++){
                for (int j=0;j<diis_size;j++){
                    B[i*d_size + j] = diis_error[i] * diis_error[j];
            }
                B[i*d_size + d_size-1] = -1;
                B[(d_size-1)*d_size+i] = -1;
                c[i] = 0.0;

            }
            B[d_size*d_size-1] = 0.0;



            c[d_size-1] = -1;

            int* IPIV = new int[d_size];
            int info,nrhs;
            nrhs =1;
            double* work;
            int lwork = -1;
            double wkopt;
            char mode[1] ={'N'};
            //dgesv_(&d_size,&nrhs,B,&d_size,IPIV,c,&d_size,&info);
            dgels_(mode,&d_size,&d_size,&nrhs,B,&d_size,c,&d_size,&wkopt,&lwork,&info);
            lwork = (int)wkopt;
            work = (double*)malloc( lwork*sizeof(double) );
            dgels_(mode,&d_size,&d_size,&nrhs,B,&d_size,c,&d_size,work,&lwork,&info);

            std::cout<<"\nC: [";
            for(int i=0;i<d_size;i++)
                std::cout<<c[i]<<", ";
            std::cout<<"]\nd_size:"<<d_size-1<<"\n\n";

            double nsum = 0.0;
            /*
            memset(Fmat,0,nroao*nroao*sizeof(double));
            for (int j=0;j<diis_size;j++){
                for(int i=0;i<nroao*nroao;i++){
                    Fmat[i]+=c[j]*diis_mat[j][i];
               }
            }*/

        }


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
