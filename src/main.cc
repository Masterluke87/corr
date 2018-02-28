#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "MATH/ops_mat.h"
#include "IO/ops_io.h"
#include "SCF/ops_rhf.h"
#include "SCF/scf.h"
#include "POSTHF/motrans.h"
#include "POSTHF/loc_mp2.h"
#include "POSTHF/can_mp2.h"
#include "POSTHF/ri.h"
#include <iomanip>
#include <omp.h>
#include <libint2.hpp>
#include <cblas.h>
#include <lapacke.h>


#define PWIDTH_L 12
#define PWIDTH_R 16



double get_integral(double* ints,long long int &istep,long long int &jstep,long long int &kstep,int  &i,int  &j,int &k,int &l){
    return ints[(i/2)*istep + (k/2)*jstep + (j/2)*kstep + l/2]*(i%2==k%2)*(j%2==l%2) -
           ints[(i/2)*istep + (l/2)*jstep + (j/2)*kstep + k/2]*(i%2==l%2)*(j%2==k%2);
}

int main(int argc, char const *argv[]) {
    int worldsize = 1;
    int rank = 0;

    /*
     *
       Part 1 : Initialization
     *
     */

    // System Information:

    systeminfo *sysinfo = new systeminfo;

    std::map<std::string,std::pair<double,double> > timer;

    timer["total"] = std::make_pair(0.0,0.0);
    timer["total"].first = omp_get_wtime();


    if(argc != 2) {
        std::cerr << "Need PREFIX\n";
        exit(1);
    }

    sysinfo->prefix = argv[1];
    read_system(sysinfo);

    sysinfo->ion_rep =  calc_ion_rep(sysinfo);
    print_header(sysinfo);


    libint2::initialize();

    libint2::BasisSet obs(sysinfo->basisNameOB,sysinfo->atoms);
    libint2::BasisSet dfbs(sysinfo->basisNameRI,sysinfo->atoms);
    obs.set_pure(true);dfbs.set_pure(true);

    /*
     *
     * PART2 : One Electron & Electron Matrices and Integrals
     *
     */

    OEints* onemats = new OEints;
    TEints* twomats = new TEints;

    allocate_onemats(sysinfo,onemats);
    // Always read the One-Electron Matrices
    // Maybe on a day in a far far future, we calculate the OEI

    read_oei(sysinfo,onemats);


    //Two possibilities:
    //1) ERI provided -> read them
    //2) NO ERI provided -> Transform Hmat to Libint -> Calculate -> ERI

    if (sysinfo->nrofint > 0) {
        std::cout << "Two electron integrals provided, reading in .." << '\n';
        twomats->intval        = new double[sysinfo->nrofint];                   //two electron integrals
        twomats->intnums       = new unsigned short[sysinfo->nrofint*4];         //two electron indices

        //ints_start = omp_get_wtime();
        read_tei(sysinfo,twomats);
        //ints_end = omp_get_wtime();
    }else{
        std::cout << "No Tei's provided ... using libint to calculate them" << '\n';

        calculate_libint_oei(sysinfo,obs,onemats);
        //ints_start = omp_get_wtime();
        calculate_libint_tei(sysinfo,obs,twomats);
        //ints_end = omp_get_wtime();

    }
    calc_S12(sysinfo,onemats);

    /*
     *
     * Part 3: Read Guess & Perform SCF
     *
     */

    /*
    read_wav_HF(prefix+".ahfw",nroao,MOens,MOs);
    {             //check if gd orbitals are provided;
        double modiag=0.0;
        for (size_t i = 0; i < nroao; i++) {
            modiag += MOs[i*nroao+i];
        }
        if (std::fabs(modiag) < 0.1) {
            std::cout << "No converged MOs provided... doing core guess" << '\n';
            for (int i =0; i<nroao*nroao; i++)
                Fmat[i] = Hmat[i];
            double *tmpmat = new double[nroao*nroao];
            diag_Fmat(nroao, Fmat,MOs,MOens,Som12, tmpmat);
            delete[] tmpmat;
        }
    }
    timer["SCF"] = std::make_pair(0.0,0.0);
    timer["SCF"].first =omp_get_wtime();
    run_scf(nroao,nroe,MOs,Pmat,Hmat,Smat,Fmat,MOens,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);
    timer["SCF"].second=omp_get_wtime();


    /*
     * PART 4:
     * Post Hartree-Fock Part:
     * BPQ of /c basis
     * transform to Bia
     * perform mp2
     *
     */
/*
    double* FMo  = new double[nroao*nroao];
    build_FMo(nroao,Fmat,FMo,MOs);
    double* prec_ints;
    {

        double trafo_start = omp_get_wtime();
        MOtrans(MOs,nroao,nroe,nrofint,sortcount,intval,intnums,&prec_ints);
        double trafo_end = omp_get_wtime();

        double mp2_start = omp_get_wtime();
        run_canonical_mp2(nroe,nroao,prec_ints,FMo);
        double mp2_end   = omp_get_wtime();
    }



    //RI-MP2 //trafo

    double* BPQ = new double[naux_2*nroao*nroao];
    calculate_ri(obs,dfbs,BPQ);
    double ritrafo_start = 0.0;
    double ritrafo_end   = 0.0;
    double rimp2_start   = 0.0;
    double rimp2_end     = 0.0;
    {
        ritrafo_start = omp_get_wtime();
        double* Bia = new double[naux_2*nroe/2*(nroao-nroe/2)];
        //read_transform_ri(prefix,nroe,nroao,naux_2,MOs,Bia);
        transform_ri(nroe,nroao,naux_2,MOs,BPQ,Bia);
        ritrafo_end = omp_get_wtime();
        rimp2_start = omp_get_wtime();
        run_canonical_mp2_ri(nroe,nroao,naux_2,Bia,FMo);
        rimp2_end = omp_get_wtime();
        delete[] Bia;
    } //end RI-MP2



    /* What the fuck!! do we need ?
     * FOCK Matrix in spin-orbital basis
     * t1 - amplides
     * t2 - amplitudes
     * eri tensor in spin orbital basis
     */

/*
    std::cout<<"nel:"<<nroe<<std::endl;
    std::cout<<"nocc:"<<nroe<<std::endl;
    std::cout<<"nvir:"<<(2*nroao - nroe)<<std::endl;
    std::cout<<"nroao:"<<nroao<<std::endl;
    std::cout<<"nob:"<<2*nroao<<std::endl;
    std::cout<<"nalpha:"<<nroao<<std::endl;
    std::cout<<"nbeta:"<<nroao<<std::endl;


    int nocc = nroe;
    int nvir = (2*nroao - nroe);
    int norb = 2*nroao;



    double* f    = new double[norb*norb];
    double* h_so = new double[norb*norb];
    memset(h_so,0,norb*norb*sizeof(double));
    memset(f,0,norb*norb*sizeof(double));

    double twoe = 0.0;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    for(int i=0; i<norb; i++)
        for(int j=0; j<norb; j++)
            for(int mu=0; mu<nroao; mu++)
                for(int nu=0; nu<nroao; nu++) {
                    h_so[i*norb+j] += MOs[(i/2)*nroao +mu]*Hmat[mu*nroao + nu]*MOs[(j/2)*nroao +nu]*(j%2==i%2);
                }



    for(int i=0; i<norb; i++)
        for(int a=0; a<norb; a++)
        {
            twoe=0.0;
            for(int j=0; j<nocc; j++)
                twoe +=get_integral(prec_ints,istep,jstep,kstep, i, j , a , j ) ;
            f[i*norb+a] = h_so[i*norb+a] + twoe;
        }
    double Eel = 0.0;
    twoe = 0.0;
    for(int i=0; i<nocc; i++) {
        for(int j=0; j<nocc; j++)
            twoe += prec_ints[(i/2)*istep + (i/2)*jstep + (j/2)*kstep + (j/2)] - prec_ints[(i/2)*istep + (j/2)*jstep + (i/2)*kstep + (j/2)]*(i%2==j%2);
        Eel += h_so[i*norb+i];
    }
    Eel+=0.5*twoe;
    double onee = 0.0;
    double fsum = 0.0;
    for(int i=0; i<nocc; i++) {
        onee += h_so[i*norb+i];
        fsum += f[i*norb+i];
    }
    std::cout<<"El:"<<Eel+ion_rep<<std::endl;
    std::cout<<"one:"<<onee+ion_rep<<std::endl;
    std::cout<<"fsum-0.5twoe:"<<fsum-(0.5*twoe)+ion_rep;



    std::cout<<"\nGenerating Guess!";

    double* T2    = new double[nocc*nocc*nvir*nvir];
    double* T2n   = new double[nocc*nocc*nvir*nvir];
    double* T1    = new double[nocc*nvir];
    double* T1n   = new double[nocc*nvir];

    memset(T1,0,nocc*nvir*sizeof(double));
    memset(T1n,0,nocc*nvir*sizeof(double));
    memset(T2n,0,nocc*nocc*nvir*nvir*sizeof(double));

    int Ti,Tj,Ta,Tb = 0;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;
    int tmpA = 0;
    int tmpB = 0;

    for (Ti = 0; Ti < nocc; Ti++) {
        for (Tj = 0; Tj < nocc; Tj++) {
            for (Ta = 0; Ta < nvir; Ta++) {
                for (Tb = 0; Tb < nvir; Tb++) {
                    tmpA = Ta+nocc;
                    tmpB = Tb+nocc;
                    T2[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] = get_integral(prec_ints,istep,jstep,kstep, Ti, Tj , tmpA , tmpB )/
                                                                    (f[(Ti)*norb + Ti] + f[Tj*norb + Tj] - f[(Ta+nocc)*norb + (Ta+nocc)] - f[(Tb+nocc)*norb + (Tb+nocc)]);

                }
            }
        }
    }

    double EMP2 = 0.0;
    for (Ti = 0; Ti < nocc; Ti++) {
        for (Tj = 0; Tj < nocc; Tj++) {
            for (Ta = 0; Ta < nvir; Ta++) {
                for (Tb = 0; Tb < nvir; Tb++) {
                    tmpA=Ta+nocc;
                    tmpB=Tb+nocc;
                    EMP2 += get_integral(prec_ints,istep,jstep,kstep, Ti, Tj , tmpA , tmpB ) *  T2[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb];
                }
            }
        }
    }

    std::cout<<"\nEMP2:"<<0.25*EMP2;






    //Intermediates:
    double* Fae    = new double [nvir*nvir];
    double* Fmi    = new double [nocc*nocc];
    double* Fme    = new double [nocc*nvir];

    double* Wmnij  = new double [nocc*nocc*nocc*nocc];
    double* Wabef  = new double [nvir*nvir*nvir*nvir];
    double* Wmbej  = new double [nocc*nvir*nvir*nocc];

    double* tau_s  = new double [nocc*nocc*nvir*nvir];
    double* tau    = new double [nocc*nocc*nvir*nvir];


    memset(Fae,0,nvir*nvir*sizeof(double));
    memset(Fmi,0,nocc*nocc*sizeof(double));
    memset(Fme,0,nocc*nvir*sizeof(double));


    double Eold = 0.0;
    double Ecc = 0.0;
    double dE   = 0.0;

    double sum = 0.0;
    int tmpc = 0;
    int tmpd = 0;
    int tmpa = 0;
    int tmpb = 0;
    int tmpe = 0;
    int tmpf = 0;

    int m_step,n_step,a_step,b_step,e_step,i_step;

    for (int iter=0;iter<20;iter++)
    {


    //f^{k}_{c} t^{c}_{k}
    sum = 0.0;
    Eold = Ecc;
    Ecc = 0.0;
    for(int k=0;k < nocc; k++)
      for(int c=0;c < nvir; c++)
        {
        tmpc=c+nocc;
        sum+= f[ k * norb + tmpc ]*T1[ k * nvir + c ];
        }
    sum*= 1.0;
    Ecc+=sum;

    //\frac{t^{dc}_{kl} v^{kl}_{dc}}{4}
    sum = 0.0;
    for(int d=0;d < nvir; d++)
      for(int c=0;c < nvir; c++)
        for(int k=0;k < nocc; k++)
          for(int l=0;l < nocc; l++)
            {
            tmpd=d+nocc;
            tmpc=c+nocc;
            sum+= T2[ k *Ti_step + l *Tj_step + d *Ta_step + c ]*get_integral(prec_ints,istep,jstep,kstep, k , l , tmpd , tmpc );
            }
    sum*= 0.250000000000000;
    Ecc+=sum;

    //\frac{t^{c}_{l} t^{d}_{k}}{2} v^{kl}_{dc}
    sum = 0.0;
    for(int c=0;c < nvir; c++)
      for(int l=0;l < nocc; l++)
        for(int d=0;d < nvir; d++)
          for(int k=0;k < nocc; k++)
            {
            tmpc=c+nocc;
            tmpd=d+nocc;
            sum+= T1[ l * nvir + c ]*T1[ k * nvir + d ]*get_integral(prec_ints,istep,jstep,kstep, k , l , tmpd , tmpc );
            }
    sum*= 0.500000000000000;
    Ecc+=sum;

    dE = Eold-Ecc;
    std::cout<<iter <<"E(CCSD)= "<<Ecc<<" dE:"<<dE<<"  "<<"\n";

    //Tau intermediates
    for (int i=0;i<nocc;i++)
        for (int j=0;j<nocc;j++)
            for (int a=0;a<nvir;a++)
                for (int b=0;b<nvir;b++){
                    tau_s[i*Ti_step + j*Tj_step + a*Ta_step + b] = T2[i*Ti_step + j*Tj_step + a*Ta_step + b] +
                                                                0.5*(T1[i*nvir + a]*
                                                                     T1[j*nvir + b] -
                                                                     T1[i*nvir + b]*
                                                                     T1[j*nvir + a] );
                    tau[i*Ti_step + j*Tj_step + a*Ta_step + b] = T2[i*Ti_step + j*Tj_step + a*Ta_step + b] +
                                                                 T1[i*nvir + a]*T1[j*nvir + b] -
                                                                 T1[i*nvir + b]*T1[j*nvir + a];

        }

    //Fae
    for (int a=0;a<nvir;a++)
        for(int e=0;e<nvir;e++){
            Fae[a*nvir +e] = (1-(a==e))*f[(a+nocc)*norb + (nocc+e)];
            for(int m=0;m<nocc;m++){
                Fae[a*nvir +e] -= 0.5* f[m*norb +(nocc+a)]*T1[m*nvir+a];
            }
            for(int m=0;m<nocc;m++)
                for(int f=0;f<nvir;f++){
                    tmpa = a+nocc;
                    tmpe = e+nocc;
                    tmpf = f+nocc;
                    Fae[a*nvir +e] += T1[m*nvir+f]*get_integral(prec_ints,istep,jstep,kstep, m , tmpa , tmpf , tmpe );
                }
            for(int m=0;m<nocc;m++)
                for(int n=0;n<nocc;n++)
                    for(int f=0;f<nvir;f++){
                         tmpe = e+nocc;
                         tmpf = f+nocc;
                         Fae[a*nvir +e] -= 0.5 * tau_s[m*Ti_step + n*Tj_step + a*Ta_step + f] *
                                            get_integral(prec_ints,istep,jstep,kstep, m , n , tmpe , tmpf );
            }
        }

    //Fmi
    for (int m=0;m<nocc;m++)
        for(int i=0;i<nocc;i++){
            Fmi[m*nocc + i] = (1-(m==i))*f[m*norb+i];
            for(int e=0;e<nvir;e++){
                Fmi[m*nocc + i] += 0.5 * T1[i*nvir + e] * f[m*norb + (nocc+e)];
            }
            for (int n=0;n<nocc;n++)
                for(int e=0;e<nvir;e++){
                    tmpe = e+nocc;
                    Fmi[m*nocc + i] += T1[n*nvir + e]*get_integral(prec_ints,istep,jstep,kstep, m , n , i , tmpe);
                }
            for (int n=0;n<nocc;n++)
                for (int e=0;e<nvir;e++)
                    for (int f=0;f<nvir;f++){
                        tmpe = e+nocc;
                        tmpf = f+nocc;
                        Fmi[m*nocc + i] += 0.5 * tau_s[i*Ti_step+n*Tj_step+e*Ta_step +f]*get_integral(prec_ints,istep,jstep,kstep, m , n ,tmpe,tmpf);
                    }

        }


    //Fme
    for(int m=0;m<nocc;m++)
        for(int e=0;e<nvir;e++){
            Fme[m*nvir + e] = f[m*norb + (nocc+e)];
            for(int n=0;n<nocc;n++)
                for(int f=0;f<nvir;f++){
                    tmpe = e+nocc;
                    tmpf = f+nocc;
                    Fme[m*nvir+e] += T1[n*nvir+f] *get_integral(prec_ints,istep,jstep,kstep, m , n , tmpe , tmpf);
                }
        }


    //Wmnij
    m_step = nocc*nocc*nocc;
    n_step = nocc*nocc;
    i_step = nocc;

    for(int m=0;m<nocc;m++)
        for(int n=0;n<nocc;n++)
            for(int i=0;i<nocc;i++)
                for(int j=0;j<nocc;j++){
                    Wmnij[m*m_step+n*n_step+i*i_step+j] = get_integral(prec_ints,istep,jstep,kstep, m , n , i , j);
                    for(int e=0;e<nvir;e++){
                        tmpe = e+nocc;
                        Wmnij[m*m_step+n*n_step+i*i_step+j] += (T1[j*nvir + e]*get_integral(prec_ints,istep,jstep,kstep, m , n , i , tmpe)-
                                                                T1[i*nvir + e]*get_integral(prec_ints,istep,jstep,kstep, m , n , j , tmpe));
                    }
                    for(int e=0;e<nvir;e++)
                        for(int f=0;f<nvir;f++)
                        {
                            tmpe = e+nocc;
                            tmpf = f+nocc;
                            Wmnij[m*m_step+n*n_step+i*i_step+j] += 0.25*tau[i*Ti_step+j*Tj_step + e*Ta_step + f]*get_integral(prec_ints,istep,jstep,kstep, m , n , tmpe , tmpf);
                        }

                }

    //Wabef

    a_step = nvir*nvir*nvir;
    b_step = nvir*nvir;
    e_step = nvir;
    for(int a=0;a<nvir;a++)
        for(int b=0;b<nvir;b++)
            for(int e=0;e<nvir;e++)
                for(int f=0;f<nvir;f++)
                {
                    tmpa = a+nocc;
                    tmpb = b+nocc;
                    tmpe = e+nocc;
                    tmpf = f+nocc;
                    Wabef[a*a_step + b*b_step + e*e_step +f] = get_integral(prec_ints,istep,jstep,kstep, tmpa , tmpb , tmpe , tmpf);

                    for (int m=0;m<nocc;m++){
                        tmpa = a+nocc;
                        tmpe = e+nocc;
                        tmpf = f+nocc;
                        Wabef[a*a_step + b*b_step + e*e_step +f] -= (T1[m*nvir+b]*get_integral(prec_ints,istep,jstep,kstep, tmpa , m , tmpe , tmpf) -
                                                                     T1[m*nvir+a]*get_integral(prec_ints,istep,jstep,kstep, tmpb , m , tmpe , tmpf));
                    }

                    for(int m=0;m<nocc;m++)
                        for(int n=0;n<nocc;n++)
                        {
                            tmpe = e+nocc;
                            tmpf = f+nocc;
                            Wabef[a*a_step + b*b_step + e*e_step +f] += 0.25*tau[m*Ti_step+n*Tj_step + a*Ta_step + b]*get_integral(prec_ints,istep,jstep,kstep, m , n , tmpe , tmpf);
                        }

                }


    //Wmbej


    m_step = nvir*nvir*nocc;
    b_step = nvir*nocc;
    e_step = nocc;

    for (int m=0;m<nocc;m++)
        for(int b=0;b<nvir;b++)
            for(int e=0;e<nvir; e++)
                for(int j=0;j<nocc;j++)
                {
                    tmpb = nocc + b;
                    tmpe = nocc + e;
                    Wmbej[m*m_step+b*b_step+e*e_step+j] = get_integral(prec_ints,istep,jstep,kstep, m , tmpb , tmpe , j);

                    for(int f=0;f<nvir; f++){
                        tmpb = b+nocc;
                        tmpe = e+nocc;
                        tmpf = f+nocc;
                        Wmbej[m*m_step+b*b_step+e*e_step+j] += T1[j*nvir+f]*get_integral(prec_ints,istep,jstep,kstep, m , tmpb , tmpe , tmpf);
                    }

                    for(int n=0;n<nocc; n++){
                        tmpe = e+nocc;
                        Wmbej[m*m_step+b*b_step+e*e_step+j] -= T1[n*nvir+b]*get_integral(prec_ints,istep,jstep,kstep, m , n , tmpe , j);
                    }
                    for(int n=0;n<nocc; n++)
                        for(int f=0;f<nvir; f++){
                            tmpe = e+nocc;
                            tmpf = f+nocc;
                            Wmbej[m*m_step+b*b_step+e*e_step+j] -= (0.5*T2[j*Ti_step+n*Tj_step+f*Ta_step+b]+T1[j*nvir+f]*T1[n*nvir+b])*get_integral(prec_ints,istep,jstep,kstep, m , n , tmpe , tmpf);
                        }

                }

    //T1n update

    for(int i=0;i<nocc;i++)
        for(int a=0;a<nvir;a++)
        {
            T1n[i*nvir+a] = f[i*norb + (nocc+a)];
            for(int e = 0;e<nvir;e++)
            {
                T1n[i*nvir+a] +=  T1[i*nvir+e] * Fae[a*nvir + e];
            }
            for(int m = 0;m<nocc;m++){
                T1n[i*nvir+a] -=  T1[m*nvir+a]*Fmi[m*nocc+i];
            }
            for(int m = 0;m<nocc;m++)
                for(int e = 0;e<nvir;e++){
                    T1n[i*nvir+a] += T2[i*Ti_step + m*Tj_step+a*Ta_step +e] *Fme[m*nvir+e];
                }
           for(int n = 0;n<nocc;n++)
               for(int f = 0;f<nvir;f++){
                   tmpa = a+nocc;
                   tmpf = f+nocc;
                   T1n[i*nvir+a] -=T1[n*nvir+f] * get_integral(prec_ints,istep,jstep,kstep, n , tmpa , i , tmpf);
               }
           for(int m = 0;m<nocc;m++)
               for(int e = 0;e<nvir;e++)
                   for(int f = 0;f<nvir;f++)
                   {
                       tmpa = a+nocc;
                       tmpf = f+nocc;
                       tmpe = e+nocc;
                       T1n[i*nvir+a] -= 0.5*(T2[i*Ti_step+m*Tj_step+e*Ta_step+f] * get_integral(prec_ints,istep,jstep,kstep, m , tmpa , tmpe , tmpf));
                   }
            for(int m = 0;m<nocc;m++)
                for(int e = 0;e<nvir;e++)
                    for(int n = 0;n<nocc;n++){
                        tmpe = e+nocc;
                        T1n[i*nvir+a] -= 0.5*(T2[m*Ti_step+n*Tj_step+a*Ta_step+e]*get_integral(prec_ints,istep,jstep,kstep, n , m , tmpe , i));
                    }

        T1n[i*nvir+a] /=( f[i*norb+i] -f[(nocc+a)*norb +nocc+a]);
        }


    //T2 Update
    double s0 = 0.0;
    double s1 = 0.0;
    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int a=0;a<nvir;a++)
                for(int b=0;b<nvir;b++){
                    tmpa = a + nocc;
                    tmpb = b + nocc;
                    T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] = get_integral(prec_ints,istep,jstep,kstep, i , j , tmpa , tmpb);
                    for(int e=0;e<nvir;e++)
                    {
                        s0 = 0.0;
                        s1 = 0.0;
                        for(int m=0;m<nocc;m++){
                            s0 += -0.5*(T1[m*nvir+b] * Fme[m*nvir+e]);
                            s1 += -0.5*(T1[m*nvir+a] * Fme[m*nvir+e]);
                        }
                        s0 += Fae[b*nvir+e];
                        s1 += Fae[a*nvir+e];
                        T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] += T2[i*Ti_step+j*Tj_step+a*Ta_step+e]*s0 - T2[i*Ti_step+j*Tj_step+b*Ta_step+e]*s1;
                    }

                    for(int m=0;m<nocc;m++)
                    {
                        s0 = 0.0;
                        s1 = 0.0;
                        for(int e=0;e<nvir;e++){
                            s0 += 0.5*(T1[j*nvir+e]*Fme[m*nvir+e]);
                            s1 += 0.5*(T1[i*nvir+e]*Fme[m*nvir+e]);
                        }
                        s0 += Fmi[m*nocc+j];
                        s1 += Fmi[m*nocc+i];
                        T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] -= (T2[i*Ti_step+m*Tj_step+a*Ta_step+b]*s0 - T2[j*Ti_step+m*Tj_step+a*Ta_step+b]*s1);

                    }
                    m_step = nocc*nocc*nocc;
                    n_step = nocc*nocc;
                    i_step = nocc;
                    for(int m=0;m<nocc;m++)
                        for(int n=0;n<nocc;n++){
                            T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] += 0.5*tau[m*Ti_step+n*Tj_step+a*Ta_step+b]*Wmnij[m*m_step+n*n_step+i*i_step+j];
                        }
                    a_step = nvir*nvir*nvir;
                    b_step = nvir*nvir;
                    e_step = nvir;
                    for(int e=0;e<nvir;e++)
                        for(int f=0;f<nvir;f++){
                            T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] += 0.5*tau[i*Ti_step+j*Tj_step+e*Ta_step+f]*Wabef[a*a_step+b*b_step+e*e_step+f];
                        }

                    m_step = nvir*nvir*nocc;
                    b_step = nvir*nocc;
                    e_step = nocc;
                    for(int m=0;m<nocc;m++)
                        for(int e=0;e<nvir;e++){
                            tmpb = b+nocc;
                            tmpa = a+nocc;
                            tmpe = e+nocc;
                            T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] +=((T2[i*Ti_step+m*Tj_step+a*Ta_step+e] * Wmbej[m*m_step+b*b_step+e*e_step+j]-T1[i*nvir+e]*T1[m*nvir+a]* get_integral(prec_ints,istep,jstep,kstep, m , tmpb , tmpe , j))-
                                                                    (T2[i*Ti_step+m*Tj_step+b*Ta_step+e] * Wmbej[m*m_step+a*b_step+e*e_step+j]-T1[i*nvir+e]*T1[m*nvir+b]* get_integral(prec_ints,istep,jstep,kstep, m , tmpa , tmpe , j)))-
                                                                   ((T2[j*Ti_step+m*Tj_step+a*Ta_step+e] * Wmbej[m*m_step+b*b_step+e*e_step+i]-T1[j*nvir+e]*T1[m*nvir+a]* get_integral(prec_ints,istep,jstep,kstep, m , tmpb , tmpe , i))-
                                                                    (T2[j*Ti_step+m*Tj_step+b*Ta_step+e] * Wmbej[m*m_step+a*b_step+e*e_step+i]-T1[j*nvir+e]*T1[m*nvir+b]* get_integral(prec_ints,istep,jstep,kstep, m , tmpa , tmpe , i)));
                     }
                    for(int e=0;e<nvir;e++){
                        tmpb = b+nocc;
                        tmpa = a+nocc;
                        tmpe = e+nocc;
                        T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] += (T1[i*nvir+e] * get_integral(prec_ints,istep,jstep,kstep, tmpa , tmpb , tmpe , j)-
                                                                 T1[j*nvir+e] * get_integral(prec_ints,istep,jstep,kstep, tmpa , tmpb , tmpe , i));
                    }
                    for(int m=0;m<nocc;m++){
                        tmpb = b+nocc;
                        tmpa = a+nocc;
                        T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] -= (T1[m*nvir+a]*get_integral(prec_ints,istep,jstep,kstep, m , tmpb , i , j))-
                                                                (T1[m*nvir+b]*get_integral(prec_ints,istep,jstep,kstep, m , tmpa , i , j) );
                    }

                T2n[i*Ti_step+j*Tj_step+a*Ta_step+b]/=(f[i*norb+i]+f[j*norb+j]-f[(nocc+a)*norb+(nocc+a)]-f[(nocc+b)*norb + nocc+b]);

                }
    /*
    //debug section
    for(int i=0;i<(nocc*nvir);i++)
        std::cout<<i<<" "<<T1n[i]<<std::endl;


    std::cout<<"ITER:"<<iter;
    std::cout<<"Fae\n";
    for(int i=0;i<nvir;i++)
        for(int j=0;j<nvir;j++)
        {
            if (abs(Fae[i*nvir+j])>1E-6){
                std::cout<<i<<" "<<j<<" : "<< Fae[i*nvir+j]<<std::endl;
            }
        }
    std::cout<<"Fmi\n";
    for(int m=0;m<nocc;m++)
        for(int i=0;i<nocc;i++)
        {
            if (abs(Fmi[m*nocc+i])>1E-6){
                std::cout<<m<<" "<<i<<" : "<< Fmi[m*nocc+i]<<std::endl;
            }
        }
    std::cout<<"Fme\n";
    for(int m=0;m<nocc;m++)
        for(int e=0;e<nvir;e++)
        {
            if (abs(Fme[m*nvir+e])>1E-6){
                std::cout<<m<<" "<<e<<" : "<< Fme[m*nvir+e]<<std::endl;
            }
        }

    std::cout<<"Wmnij\n";
    m_step = nocc*nocc*nocc;
    n_step = nocc*nocc;
    i_step = nocc;

    for(int m=0;m<nocc;m++)
        for(int n=0;n<nocc;n++)
            for(int i=0;i<nocc;i++)
                for(int j=0;j<nocc;j++){
                    if (abs(Wmnij[m*m_step+n*n_step+i*i_step+j])>1E-6){
                        std::cout<<m<<" "<<n<<" "<<i<<" "<<j<<" : "<< Wmnij[m*m_step+n*n_step+i*i_step+j]<<std::endl;
                    }
                }
    std::cout<<"Wabef\n";
    a_step = nvir*nvir*nvir;
    b_step = nvir*nvir;
    e_step = nvir;
    for(int a=0;a<nvir;a++)
        for(int b=0;b<nvir;b++)
            for(int e=0;e<nvir;e++)
                for(int f=0;f<nvir;f++){
                    if (abs(Wabef[a*a_step+b*b_step+e*e_step+f])>1E-6){
                        std::cout<<a<<" "<<b<<" "<<e<<" "<<f<<" : "<<Wabef[a*a_step+b*b_step+e*e_step+f]<<std::endl;
                    }
                }
    int wcounter =0;
    std::cout<<"Wmbej\n";
    m_step = nvir*nvir*nocc;
    b_step = nvir*nocc;
    e_step = nocc;

    for (int m=0;m<nocc;m++)
        for(int b=0;b<nvir;b++)
            for(int e=0;e<nvir; e++)
                for(int j=0;j<nocc;j++){
                    if (abs(Wmbej[m*m_step+b*b_step+e*e_step+j])>1E-6){
                        std::cout<<m<<" "<<b<<" "<<e<<" "<<j<<" : "<<Wmbej[m*m_step+b*b_step+e*e_step+j]<<std::endl;
                        wcounter++;
                    }
                }
    std::cout<<wcounter;


    std::cout<<"T1\n";
    for (int i=0;i<nocc;i++)
        for(int a=0;a<nvir;a++){
            if (abs(T1n[i*nvir+a])>1E-12){
                std::cout<<i<<" "<<a<<" : "<<T1n[i*nvir+a]<<std::endl;
            }
        }

    std::cout<<"T2\n";
    int t2counter=0;
    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int a=0;a<nvir;a++)
                for(int b=0;b<nvir;b++){
                    if (abs(T2n[i*Ti_step + j*Tj_step + a*Ta_step + b])>1E-6){
                        std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" : "<<T2n[i*Ti_step + j*Tj_step + a*Ta_step + b]<<std::endl;
                        t2counter++;
                    }
                }
    std::cout<<t2counter<<"\n";


    std::cout<<"tau\n";
    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int a=0;a<nvir;a++)
                for(int b=0;b<nvir;b++){
                    if (abs(tau[i*Ti_step + j*Tj_step + a*Ta_step + b])>1E-6){
                        std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" : "<<tau[i*Ti_step + j*Tj_step + a*Ta_step + b]<<std::endl;
                    }
                }

    std::cout<<"taus\n";
    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int a=0;a<nvir;a++)
                for(int b=0;b<nvir;b++){
                    if (abs(tau_s[i*Ti_step + j*Tj_step + a*Ta_step + b])>1E-6){
                        std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" : "<<tau_s[i*Ti_step + j*Tj_step + a*Ta_step + b]<<std::endl;
                    }
                }



    std::cout.flush();
    */
/*
    double* swapper = T2;
    T2 = T2n;
    T2n = swapper;
    swapper = T1;
    T1 = T1n;
    T1n = swapper;

    }


       /* Part 5
     * Print timings
     */



    std::cout <<"\n\nTIMINGS:\n";
    std::cout <<"--------\n";
    // std::cout << std::setw( PWIDTH_L ) << "TEI [s]:" <<std::setw( PWIDTH_R )<<ints_end - ints_start<< " s \n";
    // std::cout << std::setw( PWIDTH_L ) << "SCF [s]:" <<std::setw( PWIDTH_R )<<scf_end-scf_start<< " s \n";
    // std::cout << std::setw( PWIDTH_L ) << "4-index [s]:" <<std::setw( PWIDTH_R )<<trafo_end-trafo_start<< " s \n";
    // std::cout << std::setw( PWIDTH_L ) << "ri-4ind [s]:" <<std::setw( PWIDTH_R )<<ritrafo_end-ritrafo_start<< " s \n";
    // std::cout << std::setw( PWIDTH_L ) << "RI-MP2 [s]:" <<std::setw( PWIDTH_R )<<rimp2_end-rimp2_start<< " s \n";
    // std::cout << std::setw( PWIDTH_L ) << "MP2 [s]:" <<std::setw( PWIDTH_R )<<mp2_end-mp2_start<< " s \n";

    timer["total"].second = omp_get_wtime();
    std::cout << "------------------------------" << '\n';
    std::cout<< std::setw( PWIDTH_L ) << "Total [s]:"<<std::setw( PWIDTH_R )<<timer["total"].second - timer["total"].first<<" s \n";
    std::cout << "\n\n--END--" << '\n';
    libint2::finalize();
    return 0;
}
