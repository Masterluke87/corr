#include "ccsd.h"
#include "../UTIL/util.h"
#include "ccsd_intermediates.h"
#include <cblas.h>

double get_integral(double* ints,long long int &istep,long long int &jstep,long long int &kstep,int  &i,int  &j,int &k,int &l){
    return ints[(i/2)*istep + (k/2)*jstep + (j/2)*kstep + l/2]*(i%2==k%2)*(j%2==l%2) -
           ints[(i/2)*istep + (l/2)*jstep + (j/2)*kstep + k/2]*(i%2==l%2)*(j%2==k%2);
}

void form_fock_ri(cc_helper* CC,pHF* postHF,systeminfo* sysinfo)
{
    int norb   = CC->norb;
    int nroao  = CC->nroao;
    int nocc   = CC->nocc;
    int nvir   = CC->nvir;

    CC->f_ri      = new double[norb*norb];
    CC->h         = new double[norb*norb];

    double* f  = CC->f_ri;
    double* h  = CC->h;
    double* Hmat = CC->Hmat;
    double* MOs = CC->MOs;

    memset(h,0,norb*norb*sizeof(double));
    memset(f,0,norb*norb*sizeof(double));

    double twoe = 0.0;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    for(int i=0; i<norb; i++)
        for(int j=0; j<norb; j++)
            for(int mu=0; mu<nroao; mu++)
                for(int nu=0; nu<nroao; nu++) {
                    h[i*norb+j] += MOs[(i/2)*nroao +mu]*Hmat[mu*nroao + nu]*MOs[(j/2)*nroao +nu]*(j%2==i%2);
                }



    for(int i=0;i<nocc;i+=2){
        for(int a=0; a<nocc; a+=2){
            twoe=0.0;
            for(int j=0; j<nocc; j+=2)
                twoe += (2*cblas_ddot(sysinfo->naux_2,&(postHF->Bij[i/2*nocc/2+a/2]),nocc*nocc/4,&(postHF->Bij[j/2*nocc/2+j/2]),nocc*nocc/4)
                         -cblas_ddot(sysinfo->naux_2,&(postHF->Bij[i/2*nocc/2+j/2]),nocc*nocc/4,&(postHF->Bij[j/2*nocc/2+a/2]),nocc*nocc/4));
            f[i*norb+a] = h[i*norb+a] + twoe;
            f[(i+1)*norb+a+1] = h[i*norb+a] + twoe;
        }

    }
/*
    for(int i=0;i<nocc/2;i+=2){
        for(int a=0; a<nvir/2; a+=2){
            twoe=0.0;
            for(int j=0; j<nocc/2; j++)
                twoe += (2*cblas_ddot(sysinfo->naux_2,&(postHF->Bia[i*nvir+a]),nocc*nvir,&(postHF->Bij[j*nocc+j]),nocc*nocc)
                          -cblas_ddot(sysinfo->naux_2,&(postHF->Bij[i*nocc+j]),nocc*nocc,&(postHF->Bia[j*nvir+a]),nocc*nvir));
            f[i*norb+a+nocc] = h[i*norb+a+nocc] + twoe;
            f[(i+1)*norb+nocc+1] = h[(i+1)*norb+a+nocc+1] + twoe;
        }

    }

    for(int i=0;i<nvir/2;i+=2){
        for(int a=0; a<nocc/2; a+=2){
            twoe=0.0;
            for(int j=0; j<nocc/2; j++)
                twoe += (2*cblas_ddot(sysinfo->naux_2,&(postHF->Bia[a*nocc+i]),nocc*nvir,&(postHF->Bij[j*nocc+j]),nocc*nocc)
                          -cblas_ddot(sysinfo->naux_2,&(postHF->Bia[j*nvir+i]),nocc*nvir,&(postHF->Bij[j*nocc+a]),nocc*nocc));
            f[i*norb+a] = h[i*norb+a] + twoe;
            f[(i+1)*norb+a+1] = h[(i+1)*norb+a+1] + twoe;
        }

    }
*/
    for(int i=0;i<nvir;i+=2){
        for(int a=0; a<nvir; a+=2){
            twoe=0.0;
            for(int j=0; j<nocc; j+=2)
                twoe += (2*cblas_ddot(sysinfo->naux_2,&(postHF->Bab[a/2*nvir/2+i/2]),nvir*nvir/4,&(postHF->Bij[j/2*nocc/2+j/2]),nocc*nocc/4)
                          -cblas_ddot(sysinfo->naux_2,&(postHF->Bia[j/2*nvir/2+i/2]),nocc*nvir/4,&(postHF->Bia[j/2*nvir/2+a/2]),nocc*nvir/4));
            f[(i+nocc)*norb+nocc+a] = h[(i+nocc)*norb+nocc+a] + twoe;
            f[(i+nocc+1)*norb+a+nocc+1] = h[(i+nocc)*norb+a+nocc] + twoe;
        }

    }



/*
    for(int i=0; i<norb; i++)
        for(int a=0; a<norb; a++)
        {
            twoe=0.0;
            for(int j=0; j<nocc; j++)
                twoe +=    get_integral(postHF->prec_ints,istep,jstep,kstep, i, j , a , j ) ;
            f[i*norb+a] = h[i*norb+a] + twoe;
        }
*/

}

void form_fock_restr(cc_helper* CC,pHF* postHF){
    int nroao   = CC->nroao;
    int nocc   = CC->nocc;

    CC->f      = new double[nroao*nroao];
    CC->h      = new double[nroao*nroao];

    double* f  = CC->f;
    double* h  = CC->h;
    double* Hmat = CC->Hmat;
    double* MOs = CC->MOs;

    memset(h,0,nroao*nroao*sizeof(double));
    memset(f,0,nroao*nroao*sizeof(double));

    double twoe = 0.0;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    for(int i=0; i<nroao; i++)
        for(int j=0; j<nroao; j++)
            for(int mu=0; mu<nroao; mu++)
                for(int nu=0; nu<nroao; nu++) {
                    h[i*nroao+j] += MOs[(i)*nroao +mu]*Hmat[mu*nroao + nu]*MOs[j*nroao +nu];
                }

    for(int i=0; i<nroao; i++)
        for(int a=0; a<nroao; a++)
        {
            twoe=0.0;
            for(int j=0; j<nocc; j++)
                twoe += 2*postHF->prec_ints[i*istep + a*jstep + j*kstep +j] - postHF->prec_ints[i*istep + j*jstep + j*kstep +a];
            f[i*nroao+a] = h[i*nroao+a] + twoe;
        }
}

void form_fock(cc_helper* CC,pHF* postHF)
{
    int norb   = CC->norb;
    int nroao  = CC->nroao;
    int nocc   = CC->nocc;

    CC->f      = new double[norb*norb];
    CC->h      = new double[norb*norb];

    double* f  = CC->f;
    double* h  = CC->h;
    double* Hmat = CC->Hmat;
    double* MOs = CC->MOs;

    memset(h,0,norb*norb*sizeof(double));
    memset(f,0,norb*norb*sizeof(double));

    double twoe = 0.0;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    for(int i=0; i<norb; i++)
        for(int j=0; j<norb; j++)
            for(int mu=0; mu<nroao; mu++)
                for(int nu=0; nu<nroao; nu++) {
                    h[i*norb+j] += MOs[(i/2)*nroao +mu]*Hmat[mu*nroao + nu]*MOs[(j/2)*nroao +nu]*(j%2==i%2);
                }

    for(int i=0; i<norb; i++)
        for(int a=0; a<norb; a++)
        {
            twoe=0.0;
            for(int j=0; j<nocc; j++)
                twoe +=get_integral(postHF->prec_ints,istep,jstep,kstep, i, j , a , j ) ;
            f[i*norb+a] = h[i*norb+a] + twoe;
        }
}


void calc_e_check_ri(cc_helper* CC,pHF* postHF){
    double Eel = 0.0;
    double twoe = 0.0;
    double ion_rep = CC->ion_rep;
    int nroao = CC->nroao;
    int norb  = CC->norb;
    int nocc  = CC->nocc;

    double* h = CC->h;
    double* f = CC->f_ri;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    for(int i=0; i<nocc; i++) {
        for(int j=0; j<nocc; j++)
            twoe += postHF->prec_ints[(i/2)*istep + (i/2)*jstep + (j/2)*kstep + (j/2)] - postHF->prec_ints[(i/2)*istep + (j/2)*jstep + (i/2)*kstep + (j/2)]*(i%2==j%2);
        Eel += h[i*norb+i];
    }
    Eel+=0.5*twoe;
    double onee = 0.0;
    double fsum = 0.0;
    for(int i=0; i<nocc; i++) {
        onee += h[i*norb+i];
        fsum += f[i*norb+i];
    }
    std::cout<<"El:"<<Eel+ion_rep<<std::endl;
    std::cout<<"one:"<<onee+ion_rep<<std::endl;
    std::cout<<"fsum-0.5twoe:"<<fsum-(0.5*twoe)+ion_rep;
}


void calc_e_check_restr(cc_helper* CC,pHF* postHF){
    double Eel = 0.0;
    double twoe = 0.0;
    double ion_rep = CC->ion_rep;
    int nroao = CC->nroao;
    int nocc  = CC->nocc;

    double* prec_ints = postHF->prec_ints;
    double* h = CC->h;
    double* f = CC->f;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    for(int i=0; i<nocc; i++) {
        for(int j=0; j<nocc; j++)
            twoe += 2*prec_ints[(i)*istep + (i)*jstep + (j)*kstep + (j)] - prec_ints[(i)*istep + (j)*jstep + (i)*kstep + (j)];
        Eel += 2*h[i*nroao+i];
    }
    Eel+=twoe;
    double onee = 0.0;
    double fsum = 0.0;
    for(int i=0; i<nocc; i++) {
        onee += 2*h[i*nroao+i];
        fsum += 2*f[i*nroao+i];
    }
    std::cout<<"El:"<<Eel+ion_rep<<std::endl;
    std::cout<<"one:"<<onee+ion_rep<<std::endl;
    std::cout<<"twoe:"<<twoe<<std::endl;
    std::cout<<"fsum:"<<fsum-twoe+ion_rep;
}



void calc_e_check(cc_helper* CC,pHF* postHF){
    double Eel = 0.0;
    double twoe = 0.0;
    double ion_rep = CC->ion_rep;
    int nroao = CC->nroao;
    int norb  = CC->norb;
    int nocc  = CC->nocc;

    double* prec_ints = postHF->prec_ints;
    double* h = CC->h;
    double* f = CC->f;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    for(int i=0; i<nocc; i++) {
        for(int j=0; j<nocc; j++)
            twoe += prec_ints[(i/2)*istep + (i/2)*jstep + (j/2)*kstep + (j/2)] - prec_ints[(i/2)*istep + (j/2)*jstep + (i/2)*kstep + (j/2)]*(i%2==j%2);
        Eel += h[i*norb+i];
    }
    Eel+=0.5*twoe;
    double onee = 0.0;
    double fsum = 0.0;
    for(int i=0; i<nocc; i++) {
        onee += h[i*norb+i];
        fsum += f[i*norb+i];
    }
    std::cout<<"El:"<<Eel+ion_rep<<std::endl;
    std::cout<<"one:"<<onee+ion_rep<<std::endl;
    std::cout<<"fsum-0.5twoe:"<<fsum-(0.5*twoe)+ion_rep;
}

void ccsd_guess_restr(cc_helper *CC,pHF* postHF){

    int nvir =  CC->nvir;
    int nocc =  CC->nocc;
    int nroao = CC->nroao;

    double* prec_ints = postHF->prec_ints;
    double* f = CC->f;

    double* T1 = CC->T1;
    double* T2 = CC->T2;



    memset(T1 ,0,nocc*nvir*sizeof(double));
    memset(T2 ,0,nocc*nocc*nvir*nvir*sizeof(double));

    int Ti,Tj,Ta,Tb = 0;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;
    int tmpA = 0;
    int tmpB = 0;

    long long int nkstep = nroao;
    long long int njstep = nroao*nkstep;
    long long int nistep = nroao*njstep;


    for (Ti = 0; Ti < nocc; Ti++) {
        for (Tj = 0; Tj < nocc; Tj++) {
            for (Ta = 0; Ta < nvir; Ta++) {
                for (Tb = 0; Tb < nvir; Tb++) {
                    tmpA = Ta+nocc;
                    tmpB = Tb+nocc;
                    T2[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] = - (postHF->prec_ints[Ti*nistep + tmpA*njstep + Tj*nkstep + tmpB])/
                                                                    (-f[(Ti)*nroao + Ti] - f[Tj*nroao + Tj] + f[(Ta+nocc)*nroao + (Ta+nocc)] + f[(Tb+nocc)*nroao + (Tb+nocc)]);

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
                    EMP2 += postHF->prec_ints[Ti*nistep + tmpA*njstep + Tj*nkstep + tmpB] * (2*T2[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] - T2[Ti*Ti_step + Tj*Tj_step + Tb*Ta_step + Ta]);
                }
            }
        }
    }
    std::cout<<"\n\nEMP2: "<<EMP2;
}



void ccsd_guess(cc_helper *CC,pHF* postHF){

    int nvir =  CC->nvir;
    int nocc =  CC->nocc;
    int norb =  CC->norb;
    int nroao = CC->nroao;

    double* prec_ints = postHF->prec_ints;
    double* f = CC->f;

    double* T1 = CC->T1;
    double* T2 = CC->T2;



    memset(T1 ,0,nocc*nvir*sizeof(double));
    memset(T2 ,0,nocc*nocc*nvir*nvir*sizeof(double));

    int Ti,Tj,Ta,Tb = 0;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;
    int tmpA = 0;
    int tmpB = 0;

    long long int nkstep = norb;
    long long int njstep = norb*nkstep;
    long long int nistep = norb*njstep;


    for (Ti = 0; Ti < nocc; Ti++) {
        for (Tj = 0; Tj < nocc; Tj++) {
            for (Ta = 0; Ta < nvir; Ta++) {
                for (Tb = 0; Tb < nvir; Tb++) {
                    tmpA = Ta+nocc;
                    tmpB = Tb+nocc;
                    T2[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] = postHF->ints_so[Ti*nistep + Tj*njstep + tmpA*nkstep + tmpB]/
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
                    EMP2 +=   postHF->ints_so[Ti*nistep + Tj*njstep + tmpA*nkstep + tmpB] *  T2[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb];
                }
            }
        }
    }
    std::cout<<"\n\nEMP2: "<<0.25*EMP2;
}


void ccsd_energy_restr(cc_helper *CC, pHF* postHF){

    int nocc = CC->nocc;
    int nvir = CC->nvir;
    int nroao = CC->nroao;

    double *prec_ints = postHF->prec_ints;
    double *f = CC->f;

    double *T1 = CC->T1;
    double *T2 = CC->T2;

    CC->Ecc_old = CC->Ecc;
    CC->Ecc = 0.0;


    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    double sum = 0.0;

    int tmpa = 0;
    int tmpb = 0;
    for(int i=0; i<nocc; i++)
      for(int a=0;a<nvir; a++)
        {
        tmpa=a+nocc;
        sum+= f[ i * nroao + tmpa ]*T1[ i * nvir + a ];
        }
    sum*= 2.0;
    CC->Ecc+=sum;

    //\frac{t^{dc}_{kl} v^{kl}_{dc}}{4}
    sum = 0.0;
    for(int a=0;a < nvir; a++)
      for(int b=0;b < nvir; b++)
        for(int i=0;i < nocc; i++)
          for(int j=0;j < nocc; j++)
            {
                tmpa=a+nocc;
                tmpb=b+nocc;
                sum+= (2*postHF->prec_ints[i*istep +tmpa*jstep + j*kstep +tmpb ] - postHF->prec_ints[i*istep +tmpb*jstep + j*kstep +tmpa])
                        * (CC->T2[i*Ti_step + j*Tj_step +  a*Ta_step + b] + CC->T1[i*nvir+a]*CC->T1[j*nvir+b]);
            }

    sum*= 1.0000000000000;
    CC->Ecc+=sum;

    CC->DE = CC->Ecc_old-CC->Ecc;
    std::cout<<CC->iter <<"E(CCSD)= "<<CC->Ecc<<" dE:"<<CC->DE<<"  "<<"\n";

}


void ccsd_energy(cc_helper *CC, pHF* postHF){

    int nocc = CC->nocc;
    int nvir = CC->nvir;
    int norb = CC->norb;
    int nroao = CC->nroao;

    double *prec_ints = postHF->prec_ints;
    double *f = CC->f;

    double *T1 = CC->T1;
    double *T2 = CC->T2;

    CC->Ecc_old = CC->Ecc;
    CC->Ecc = 0.0;

    int tmpc = 0;
    int tmpd = 0;

    double* ints_so = postHF->ints_so;

    long long int kstep = norb;
    long long int jstep = norb*kstep;
    long long int istep = norb*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    double sum = 0.0;


    for(int k=0;k < nocc; k++)
      for(int c=0;c < nvir; c++)
        {
        tmpc=c+nocc;
        sum+= f[ k * norb + tmpc ]*T1[ k * nvir + c ];
        }
    sum*= 1.0;
    CC->Ecc+=sum;

    //\frac{t^{dc}_{kl} v^{kl}_{dc}}{4}
    sum = 0.0;
    for(int d=0;d < nvir; d++)
      for(int c=0;c < nvir; c++)
        for(int k=0;k < nocc; k++)
          for(int l=0;l < nocc; l++)
            {
            tmpd=d+nocc;
            tmpc=c+nocc;
            sum+= T2[ k *Ti_step + l *Tj_step + d *Ta_step + c ]*ints_so[k*istep + l*jstep + tmpd*kstep + tmpc];
            }
    sum*= 0.250000000000000;
    CC->Ecc+=sum;

    //\frac{t^{c}_{l} t^{d}_{k}}{2} v^{kl}_{dc}
    sum = 0.0;
    for(int c=0;c < nvir; c++)
      for(int l=0;l < nocc; l++)
        for(int d=0;d < nvir; d++)
          for(int k=0;k < nocc; k++)
            {
            tmpc=c+nocc;
            tmpd=d+nocc;
            sum+= T1[ l * nvir + c ]*T1[ k * nvir + d ]*ints_so[k*istep + l*jstep + tmpd*kstep + tmpc];
            }
    sum*= 0.500000000000000;
    CC->Ecc+=sum;

    CC->DE = CC->Ecc_old-CC->Ecc;
    std::cout<<CC->iter <<"E(CCSD)= "<<CC->Ecc<<" dE:"<<CC->DE<<"  "<<"\n";

}

void ccsd_build_taus(cc_helper *CC,cc_intermediates* CC_int){
    double *tau   = CC_int->tau;
    double *tau_s = CC_int->tau_s;

    double *T2 = CC->T2;
    double *T1 = CC->T1;

    int nocc = CC->nocc;
    int nvir = CC->nvir;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;


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

}


void ccsd_build_intermediates(cc_helper* CC,cc_intermediates *CC_int,pHF* postHF)
{
    build_Fae(CC,CC_int,postHF);
    build_Fmi(CC,CC_int,postHF);
    build_Fme(CC,CC_int,postHF);
    build_Wmnij(CC,CC_int,postHF);
    build_Wabef(CC,CC_int,postHF);
    build_Wmbej(CC,CC_int,postHF);
}


void ccsd_update_T1(cc_helper* CC,cc_intermediates *CC_int,pHF* postHF){
    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int nroao = CC->nroao;
    int norb  = CC->norb;

    double* T1n = CC->T1n;

    double* T1 =  CC->T1;
    double* T2 =  CC->T2;
    double*  f =  CC->f;

    double* Fae = CC_int->Fae;
    double* Fmi = CC_int->Fmi;
    double* Fme = CC_int->Fme;

    double* prec_ints = postHF->prec_ints;

    int tmpa,tmpe,tmpf;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;


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

}


void ccsd_update_T2(cc_helper* CC,cc_intermediates *CC_int,pHF* postHF){
    double s0 = 0.0;
    double s1 = 0.0;

    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int norb = CC->norb;

    double* T2n = CC->T2n;

    double* T1    = CC->T1;
    double* T2    = CC->T2;
    double* tau   = CC_int->tau;
    double*  f    = CC->f;
    double* Wmnij = CC_int->Wmnij;
    double* Wabef = CC_int->Wabef;
    double* Wmbej = CC_int->Wmbej;

    double* Fae = CC_int->Fae;
    double* Fmi = CC_int->Fmi;
    double* Fme = CC_int->Fme;

    double* ints_so = postHF->ints_so;

    long long int kstep = norb;
    long long int jstep = norb*kstep;
    long long int istep = norb*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    int m_step,n_step,i_step,a_step,e_step,b_step;
    int tmpa,tmpb,tmpe;

    int idx = 0;







    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nocc*nocc,nvir*nvir,nocc*nocc,0.5,Wmnij,nocc*nocc,tau,nvir*nvir,0.0,T2n,nvir*nvir);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,nocc*nocc,nvir*nvir,nvir*nvir,0.5,tau,nvir*nvir,Wabef,nvir*nvir,1.0,T2n,nvir*nvir);





    for(int i=0;i<nocc;i++)
        for(int j=0;j<=i;j++)
            for(int a=0;a<nvir;a++)
                for(int b=0;b<=a;b++){

                    tmpa = a + nocc;
                    tmpb = b + nocc;

                    idx = i*Ti_step+j*Tj_step+a*Ta_step+b;

                    T2n[idx] += ints_so[i*istep + j*jstep + tmpa*kstep + tmpb];

                    for(int e=0;e<nvir;e++)
                    {
                        s0 = 0.0;
                        s1 = 0.0;

                        tmpb = b+nocc;
                        tmpa = a+nocc;
                        tmpe = e+nocc;

                        T2n[idx] += (T1[i*nvir+e] * ints_so[tmpa*istep + tmpb*jstep + tmpe*kstep + j] -
                                                                 T1[j*nvir+e] * ints_so[tmpa*istep + tmpb*jstep + tmpe*kstep + i]);
                        for(int m=0;m<nocc;m++){
                            s0 += -0.5*(T1[m*nvir+b] * Fme[m*nvir+e]);
                            s1 += -0.5*(T1[m*nvir+a] * Fme[m*nvir+e]);
                        }
                        s0 += Fae[b*nvir+e];
                        s1 += Fae[a*nvir+e];
                        T2n[idx] += T2[i*Ti_step+j*Tj_step+a*Ta_step+e]*s0 - T2[i*Ti_step+j*Tj_step+b*Ta_step+e]*s1;



                        a_step = nvir*nvir*nvir;
                        b_step = nvir*nvir;
                        e_step = nvir;



                    }

                    for(int m=0;m<nocc;m++)
                    {
                        s0 = 0.0;
                        s1 = 0.0;

                        m_step = nvir*nvir*nocc;
                        b_step = nvir*nocc;
                        e_step = nocc;
                        for(int e=0;e<nvir;e++){
                            s0 += 0.5*(T1[j*nvir+e]*Fme[m*nvir+e]);
                            s1 += 0.5*(T1[i*nvir+e]*Fme[m*nvir+e]);
                            tmpb = b+nocc;
                            tmpa = a+nocc;
                            tmpe = e+nocc;
                            T2n[idx] += (T2[i*Ti_step+m*Tj_step+a*Ta_step+e] * Wmbej[m*m_step+b*b_step+e*e_step+j]
                                         -T1[i*nvir+e]*T1[m*nvir+a]*ints_so[tmpb*istep + m*jstep +j*kstep+ tmpe]
                                         -T2[i*Ti_step+m*Tj_step+b*Ta_step+e] * Wmbej[m*m_step+a*b_step+e*e_step+j]
                                         +T1[i*nvir+e]*T1[m*nvir+b]*ints_so[tmpa*istep + m*jstep +j*kstep+ tmpe]
                                         -T2[j*Ti_step+m*Tj_step+a*Ta_step+e] * Wmbej[m*m_step+b*b_step+e*e_step+i]
                                         +T1[j*nvir+e]*T1[m*nvir+a]*ints_so[tmpb*istep + m*jstep +i*kstep+ tmpe]
                                         +T2[j*Ti_step+m*Tj_step+b*Ta_step+e] * Wmbej[m*m_step+a*b_step+e*e_step+i]
                                         -T1[j*nvir+e]*T1[m*nvir+b]*ints_so[tmpa*istep + m*jstep +i*kstep+ tmpe]);

                        }
                        tmpb = b+nocc;
                        tmpa = a+nocc;
                        T2n[idx] -= (T1[m*nvir+a]*ints_so[m*istep + tmpb*jstep + i*kstep + j])-
                                                                (T1[m*nvir+b]*ints_so[m*istep + tmpa*jstep + i*kstep + j]);
                        s0 += Fmi[m*nocc+j];
                        s1 += Fmi[m*nocc+i];
                        T2n[idx] -= (T2[i*Ti_step+m*Tj_step+a*Ta_step+b]*s0 - T2[j*Ti_step+m*Tj_step+a*Ta_step+b]*s1);

                    }



                    //T2n[i*Ti_step+j*Tj_step+a*Ta_step+b] += 0.5 * cblas_ddot(nvir*nvir,&tau[i*Ti_step+j*Tj_step],1,&Wabef[a*a_step+b*b_step],1);






                T2n[i*Ti_step+j*Tj_step+a*Ta_step+b]/=(f[i*norb+i]+f[j*norb+j]-f[(nocc+a)*norb+(nocc+a)]-f[(nocc+b)*norb + nocc+b]);

                T2n[j*Ti_step+i*Tj_step+a*Ta_step+b] = -T2n[i*Ti_step+j*Tj_step+a*Ta_step+b];
                T2n[i*Ti_step+j*Tj_step+b*Ta_step+a] = -T2n[i*Ti_step+j*Tj_step+a*Ta_step+b];
                T2n[j*Ti_step+i*Tj_step+b*Ta_step+a] = T2n[i*Ti_step+j*Tj_step+a*Ta_step+b];



    }
}

void allocate_amplitudes_intermediates_restr(cc_helper* CC,cc_intermediates* CC_int){
    int nocc=CC->nocc;
    int nvir=CC->nvir;

    long int mem_a = 2* (nocc*nvir) +
                     2*(nocc*nocc*nvir*nvir);

    long int mem_int = (nvir*nvir) + (nocc*nocc) + (nocc*nvir) +
                       (nocc*nocc*nocc*nocc)+
                       (nvir*nvir*nvir*nvir)+
                       (nocc*nvir*nvir*nocc)+
                       (nocc*nocc*nvir*nvir)+
                       (nocc*nocc*nvir*nvir);

    std::cout<<"\n Allocating "<<mem_a/1024/1024 <<"MB("<<mem_a/1024/1024/1024<<"Gb) for amplitues\n";
    std::cout<<"Allocating "<<mem_int/1024/1024<<"MB("<<mem_int/1024/1024/1024<<"Gb) for intermediates\n";
    std::cout.flush();
    CC->pMem = new double[(mem_a+mem_int)];
    long int inc = 0;

    CC->T1  = &(CC->pMem[inc]);inc+=nocc*nvir;
    CC->T1n = &(CC->pMem[inc]);inc+=nocc*nvir;

    CC->T2  = &(CC->pMem[inc]);inc+=nocc*nocc*nvir*nvir;
    CC->T2n = &(CC->pMem[inc]);inc+=nocc*nocc*nvir*nvir;

   // CC_int->Fae = &(CC->pMem[inc]); inc+=nvir*nvir;
   // CC_int->Fmi = &(CC->pMem[inc]); inc+=nocc*nocc;
   // CC_int->Fme = &(CC->pMem[inc]); inc+=nocc*nvir;

  //  CC_int->Wmnij = &(CC->pMem[inc]); inc+=nocc*nocc*nocc*nocc;
  //  CC_int->Wabef = &(CC->pMem[inc]); inc+=nvir*nvir*nvir*nvir;
  //  CC_int->Wmbej = &(CC->pMem[inc]); inc+=nocc*nvir*nvir*nocc;

  //  CC_int->tau   = &(CC->pMem[inc]); inc+=nocc*nocc*nvir*nvir;
  //  CC_int->tau_s = &(CC->pMem[inc]);

}

void allocate_amplitudes_intermediates(cc_helper* CC,cc_intermediates* CC_int){
    int nocc=CC->nocc;
    int nvir=CC->nvir;

    long int mem_a = 2* (nocc*nvir) +
                     2*(nocc*nocc*nvir*nvir);

    long int mem_int = (nvir*nvir) + (nocc*nocc) + (nocc*nvir) +
                       (nocc*nocc*nocc*nocc)+
                       (nvir*nvir*nvir*nvir)+
                       (nocc*nvir*nvir*nocc)+
                       (nocc*nocc*nvir*nvir)+
                       (nocc*nocc*nvir*nvir);

    std::cout<<"\n Allocating "<<mem_a/1024/1024 <<"MB("<<mem_a/1024/1024/1024<<"Gb) for amplitues\n";
    std::cout<<"Allocating "<<mem_int/1024/1024<<"MB("<<mem_int/1024/1024/1024<<"Gb) for intermediates\n";
    std::cout.flush();
    CC->pMem = new double[(mem_a+mem_int)];
    long int inc = 0;

    CC->T1  = &(CC->pMem[inc]);inc+=nocc*nvir;
    CC->T1n = &(CC->pMem[inc]);inc+=nocc*nvir;

    CC->T2  = &(CC->pMem[inc]);inc+=nocc*nocc*nvir*nvir;
    CC->T2n = &(CC->pMem[inc]);inc+=nocc*nocc*nvir*nvir;

    CC_int->Fae = &(CC->pMem[inc]); inc+=nvir*nvir;
    CC_int->Fmi = &(CC->pMem[inc]); inc+=nocc*nocc;
    CC_int->Fme = &(CC->pMem[inc]); inc+=nocc*nvir;

    CC_int->Wmnij = &(CC->pMem[inc]); inc+=nocc*nocc*nocc*nocc;
    CC_int->Wabef = &(CC->pMem[inc]); inc+=nvir*nvir*nvir*nvir;
    CC_int->Wmbej = &(CC->pMem[inc]); inc+=nocc*nvir*nvir*nocc;

    CC_int->tau   = &(CC->pMem[inc]); inc+=nocc*nocc*nvir*nvir;
    CC_int->tau_s = &(CC->pMem[inc]);

}

void calc_ints_so(cc_helper* CC,pHF* postHF)
{
    int norb = CC->norb;
    int nroao = CC->nroao;

    postHF->ints_so = new double[norb*norb*norb*norb];



    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;


    long long int nkstep = norb;
    long long int njstep = norb*nkstep;
    long long int nistep = norb*njstep;

    for(int i=0;i<norb;i++)
        for(int j=0;j<norb;j++)
            for(int k=0;k<norb;k++)
                for(int l=0;l<norb;l++)
                    postHF->ints_so[i*nistep+j*njstep+k*nkstep+l] = get_integral(postHF->prec_ints,istep,jstep,kstep,i,j,k,l);



}


void ccsd_ur(systeminfo* sysinfo,OEints* onemats,pHF* postHF){
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
    cc_intermediates * CC_int = new cc_intermediates;

    std::cout << "\nCCSD (Spin-Orbital Formulation ):\n-----------" << '\n';


    calc_ints_so(CC,postHF);
    allocate_amplitudes_intermediates(CC,CC_int);

    form_fock(CC,postHF);

   // form_fock_ri(CC,postHF,sysinfo);

    calc_e_check(CC,postHF);
    ccsd_guess(CC,postHF);


    for (CC->iter=0;CC->iter<20;CC->iter++)
    {

        ccsd_energy(CC,postHF);
        ccsd_build_taus(CC,CC_int);
        ccsd_build_intermediates(CC,CC_int,postHF);
        ccsd_update_T1(CC,CC_int,postHF);
        ccsd_update_T2(CC,CC_int,postHF);


        double* swapper = CC->T2;
        CC->T2  = CC->T2n;
        CC->T2n = swapper;
        swapper = CC->T1;
        CC->T1  = CC->T1n;
        CC->T1n = swapper;
    }

    if (abs(CC->Ecc - (-0.2232340077))>1E-10){
        std::cout<<abs(CC->Ecc - (-0.2232340077))<<"!!!!!!!!!!!!!!!!!!!!!Falue?";
    }
    else{
        std::cout<<"GD JOB";
    }

}



void ccsd_restr(systeminfo* sysinfo,OEints* onemats,pHF* postHF){
    cc_helper *CC = new cc_helper;
    CC->nocc  = sysinfo->nroe/2;
    CC->nroao = sysinfo->nroao;
    CC->nvir  = (sysinfo->nroao - CC->nocc);
    CC->MOs   = onemats->MOs;
    CC->Hmat  = onemats->Hmat;
    CC->ion_rep = sysinfo->ion_rep;
    CC->Ecc     = 0.0;
    CC->Ecc_old = 0.0;
    CC->DE      = 100000.0;
    CC->iter    = 0;
    cc_intermediates_restr * CC_int = new cc_intermediates_restr;

    std::cout << "\nCCSD (RESTRICTED Formulation):\n------------------------------" << '\n';

    int nocc = CC->nocc;
    int nvir = CC->nvir;
    int nroao =CC->nroao;

    CC->T1 = new double[nocc*nvir];
    CC->T2 = new double[nocc*nocc*nvir*nvir];

    form_fock_restr(CC,postHF);
    for (int i=0;i<sysinfo->nroao;i++){
        std::cout<<CC->f[i*CC->nroao+i] << "   "<<CC->f[i*CC->nroao + i]<<std::endl;

    }
    //allocate_amplitudes_intermediates_restr(CC,CC_int);
    calc_e_check_restr(CC,postHF);
    ccsd_guess_restr(CC,postHF);






    CC_int->Hik    = new double [nocc*nocc];
    CC_int->Hca    = new double [nvir*nvir];
    CC_int->Hck    = new double [nocc*nvir];


    double *Hik = CC_int->Hik;
    double *Hca = CC_int->Hca;
    double *Hck = CC_int->Hck;

    double* tau = new double[nocc*nocc*nvir*nvir];
    double* Gik = new double[nocc*nocc];
    double* Gca = new double[nvir*nvir];

    double* Aijkl = new double[nocc*nocc*nocc*nocc];
    double* Bcdab = new double[nvir*nvir*nvir*nvir];
    double* Jicak = new double[nocc*nvir*nvir*nocc];
    double* Kicka = new double[nocc*nvir*nocc*nvir];
    double* r1    = new double[nocc*nvir];
    double* r2    = new double[nocc*nocc*nvir*nvir];
    double* Pijab = new double[nocc*nocc*nvir*nvir];

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int K_kstep = nvir;
    int K_cstep = K_kstep*nocc;
    int K_istep = K_cstep*nvir;

    int A_kstep = nocc;
    int A_jstep = A_kstep*nocc;
    int A_istep = A_jstep*nocc;

    int B_astep = nvir;
    int B_dstep = B_astep*nvir;
    int B_cstep = B_dstep*nvir;

    int J_astep = nocc;
    int J_cstep = J_astep*nvir;
    int J_istep = J_cstep*nvir;

    for(int iter=0;iter<10;iter++){

    ccsd_energy_restr(CC,postHF);

    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int a=0;a<nvir;a++)
                for(int b=0;b<nvir;b++){
                    //tau[ijab] = t2[ijab]+t1[ia]*t1[jb]
                    tau[i*Ti_step+j*Tj_step+a*Ta_step+b] = CC->T2[i*Ti_step+j*Tj_step+a*Ta_step+b] + CC->T1[i*nvir+a]* CC->T1[j*nvir+b];
                }

    int tmpc,tmpd;
    for(int i=0;i<nocc;i++)
        for(int k=0;k<nocc;k++){
            Hik[i*nocc+k] = 0.0;
            for(int l=0;l<nocc;l++)
                for(int c=0;c<nvir;c++)
                    for(int d=0;d<nvir;d++){
                        tmpc = c+nocc;
                        tmpd = d+nocc;
                        //Hik[ik] = [2(kcld)-(kdlc)]tau[ilcd]
                        Hik[i*nocc+k] += (2*postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd] - postHF->prec_ints[k*istep + tmpd*jstep + l*kstep + tmpc])
                                         *tau[i*Ti_step + l*Tj_step +  c*Ta_step +  d];
                    }

        }

    double sum = 0.0;
    for(int a=0;a<nvir;a++)
        for(int c=0;c<nvir;c++){
            Hca[c*nvir+a] = 0.0;
            sum = 0.0;

                for(int k=0;k<nocc;k++)
                    for(int l=0;l<nocc;l++)
                         for(int d=0;d<nvir;d++){
                             tmpc = c+nocc;
                             tmpd = d+nocc;
                             //[2*(kcld) - (kdlc)]tau[klad]
                             sum += (2*postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd]
                                     - postHF->prec_ints[k*istep + tmpd*jstep + l*kstep + tmpc])
                                       *tau[k*Ti_step + l*Tj_step +  a*Ta_step +  d];

                         }

            Hca[c*nvir+a] =-sum;
        }

    for(int c=0;c<nvir;c++)
        for(int k=0;k<nocc;k++){
            Hck[c*nocc+k] = 0.0;
                for(int l=0;l<nocc;l++)
                    for(int d=0;d<nvir;d++){
                        tmpc = c+nocc;
                        tmpd = d+nocc;
                        //H[ck] = [2*(kcld)-(kdlc)]T1[ld]
                        Hck[c*nocc+k]+= (2*postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd]
                                         - postHF->prec_ints[k*istep + tmpd*jstep + l*kstep + tmpc])*CC->T1[l*nvir+d];

                    }
        }



    for (int i=0;i<nocc;i++)
        for(int k=0;k<nocc;k++){
            Gik[i*nocc + k] = Hik[i*nocc+k];
            for(int c=0;c<nvir;c++)
                for(int l=0;l<nocc;l++){
                    tmpc = c+nocc;
                    //Gik[ik]=[2*(ikcl) - (ilck)]*T1[lc]
                    Gik[i*nocc + k] += (2*postHF->prec_ints[i*istep + k*jstep + tmpc*kstep + l]
                                        - postHF->prec_ints[i*istep + l*jstep + tmpc*kstep + k])*CC->T1[l*nvir+c];
                }
        }


    int tmpa,tmpb;
    for (int a=0;a<nvir;a++)
        for(int c=0;c<nvir;c++){
            Gca[c*nvir+a] = Hca[c*nvir+a];
                for(int d=0;d<nvir;d++)
                    for(int k=0;k<nocc;k++){
                        tmpa = a+nocc;
                        tmpc = c+nocc;
                        tmpd = d+nocc;
                        //Gca[ca] += [2(acdk)- (adck)]T1[kd]
                        Gca[c*nvir+a] += (2*postHF->prec_ints[tmpa*istep + tmpc*jstep + tmpd*kstep + k]
                                          - postHF->prec_ints[tmpa*istep + tmpd*jstep + tmpc*kstep + k])*CC->T1[k*nvir+d];
                    }
        }




    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int k=0;k<nocc;k++)
                for(int l=0;l<nocc;l++){
                    //A[ijkl] = (ikjl)
                    Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] = postHF->prec_ints[i*istep + k*jstep + j*kstep + l];
                    for(int c=0;c<nvir;c++){
                        tmpc = c+nocc;
                        //A[ijkl] += (ikcl)T[jc]
                        Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] += postHF->prec_ints[i*istep + k*jstep + tmpc*kstep + l]*CC->T1[j*nvir+c];
                    }
                    for(int c=0;c<nvir;c++){
                        tmpc = c + nocc;
                        //A[ijkl] += (ckjl)T[ic]
                        Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] += postHF->prec_ints[tmpc*istep + k*jstep + j*kstep + l]*CC->T1[i*nvir+c];
                    }
                    for(int c=0;c<nvir;c++)
                        for(int d=0;d<nvir;d++){
                            tmpc = c + nocc;
                            tmpd = d + nocc;
                            //A[ijkl]+= (kcld)*tau[ijcd]
                            Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] += postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd]*tau[i*Ti_step + j*Tj_step + c*Ta_step +d];

                    }
                }



    for(int c=0; c<nvir;c++)
        for(int d=0; d<nvir;d++)
            for(int a=0; a<nvir;a++)
                for(int b=0; b<nvir;b++){
                    tmpa = a + nocc;
                    tmpb = b + nocc;
                    tmpc = c + nocc;
                    tmpd = d + nocc;
                    //B[cdab] = (acbd)
                    Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b] = postHF->prec_ints[tmpa*istep + tmpc*jstep + tmpb*kstep + tmpd];
                    for (int k=0;k<nocc;k++){
                        //B[cdab] -= (acdk)T[kb]
                        Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b] -=  postHF->prec_ints[tmpa*istep + tmpc*jstep + tmpd*kstep + k]*CC->T1[k*nvir +b];  //(acdk)t1[kb]
                        //B[cdab] -= (bdck)T[ka]
                        Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b] -=  postHF->prec_ints[tmpb*istep + tmpd*jstep + tmpc*kstep + k]*CC->T1[k*nvir +a];  //(bdck)t1[ka]
                    }



    }



    for(int i=0;i<nocc;i++)
        for(int c=0;c<nvir;c++)
            for(int a=0;a<nvir;a++)
                for(int k=0;k<nocc;k++){
                    tmpc = c+nocc;
                    tmpa = a+nocc;
                    //J[icak] = (aick)
                    Jicak[i*J_istep + c*J_cstep + a*J_astep +k] = postHF->prec_ints[tmpa*istep + i*jstep + tmpc*kstep + k];

                    for (int l=0;l<nocc;l++){
                        //J[icak] -= (ilck)T[la]
                        Jicak[i*J_istep + c*J_cstep + a*J_astep +k] -= postHF->prec_ints[i*istep + l*jstep + tmpc*kstep + k]*
                                                                          CC->T1[l*nvir+a];
                    }

                    for (int d=0;d<nvir;d++){
                        tmpd = d + nocc;
                        //J[icak] += (adck)T[id}
                        Jicak[i*J_istep + c*J_cstep + a*J_astep +k] += postHF->prec_ints[tmpa*istep + tmpd*jstep + tmpc*kstep + k]
                                                                *CC->T1[i*nvir+d];
                    }

                    for (int d=0;d<nvir;d++)
                        for (int l=0;l<nocc;l++){
                            tmpd = d + nocc;
                            //J[icak] -= 1/2(ckdl)(t2[ilda]+2*t[id]*t[la]
                            Jicak[i*J_istep + c*J_cstep + a*J_astep +k] -= 0.5*postHF->prec_ints[tmpc*istep + k*jstep + tmpd*kstep + l]*
                                                                           (CC->T2[i*Ti_step + l*Tj_step + d*Ta_step +a] +
                                                                          2*CC->T1[i*nvir+d] * CC->T1[l*nvir+a]);
                        }

                    for (int d=0;d<nvir;d++)
                        for (int l=0;l<nocc;l++){
                            tmpd = d + nocc;
                            //J[icak] += 1/2[2(ckdl)-(dkcl)]t2[ilad]
                            Jicak[i*J_istep + c*J_cstep + a*J_astep +k] += 0.5*(2*postHF->prec_ints[tmpc*istep + k*jstep + tmpd*kstep + l]-
                                                                                  postHF->prec_ints[tmpd*istep + k*jstep + tmpc*kstep + l])*CC->T2[i*Ti_step + l*Tj_step + a*Ta_step +d];


                        }


                }



    for(int i=0;i<nocc;i++)
        for(int c=0;c<nvir;c++)
            for(int k=0;k<nocc;k++)
                for(int a=0;a<nvir;a++)
                {
                    tmpa = a + nocc;
                    tmpc = c + nocc;
                    //K[icka] = (ikac)
                    Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] = postHF->prec_ints[i*istep + k*jstep + tmpa*kstep + tmpc];
                    for(int l=0;l<nocc;l++){
                        //K[icka] -= (ikcl)T1[la]
                        Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] -= postHF->prec_ints[i*istep + k*jstep + tmpc*kstep + l]*CC->T1[l*nvir+a];
                    }
                    for(int d=0;d<nvir;d++){
                        // K[icka]+= (dkac)*T1[id]
                        tmpd = d+nocc;
                        Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] += postHF->prec_ints[tmpd*istep + k*jstep + tmpa*kstep + tmpc]*CC->T1[i*nvir+d];
                    }
                    for (int d=0;d<nvir;d++)
                        for (int l=0;l<nocc;l++){
                            tmpd = d + nocc;
                            //K[icka] -= 1/2(dkcl)(t2[ilda]*2T1[id]T1[la]
                            Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] -= 0.5*postHF->prec_ints[tmpd*istep + k*jstep + tmpc*kstep + l]*
                                                                            (CC->T2[i*Ti_step + l*Tj_step + d*Ta_step+ a]
                                                                         + 2*CC->T1[i*nvir+d]*CC->T1[l*nvir+a]);


                        }

                }


    //update amplitudes


    memset(r1   ,0 ,nocc*nvir*sizeof(double));
    memset(r2   ,0 ,nocc*nocc*nvir*nvir*sizeof(double));
    memset(Pijab,0 ,nocc*nocc*nvir*nvir*sizeof(double));

    double* f = CC->f;

    for(int i=0;i<nocc;i++)
        for(int a=0;a<nvir;a++){
            for(int c=0;c<nvir;c++){
                //r[ia] += Hca[ca]T1[ic]
                r1[i*nvir +a] += Hca[c*nvir+a]*CC->T1[i*nvir+c];
            }
            for(int k=0;k<nocc;k++){
                //r[ia] -= Hik[ik]T1[ka]
                 r1[i*nvir +a] -= Hik[i*nocc+k]*CC->T1[k*nvir+a];
            }
            for(int c=0;c<nvir;c++)
                for(int k=0;k<nocc;k++){
                    //r[ia]+= Hck[ck]*(2*T2[kica] - T2[ikca] + T1[ic]*T1[ka])
                    r1[i*nvir +a] +=  Hck[c*nocc+k]*
                                      (2*CC->T2[k*Ti_step + i*Tj_step + c*Ta_step + a]
                                       - CC->T2[i*Ti_step + k*Tj_step + c*Ta_step + a]
                                       + CC->T1[i*nvir+c]*CC->T1[k*nvir+a]);
                }
            for(int c=0;c<nvir;c++)
                for(int k=0;k<nocc;k++){
                       tmpc = c+nocc;
                       tmpa = a+nocc;
                       // r[ia]+= [2(ckai)-(ikac)]T1[kc]
                       r1[i*nvir +a] += (2*postHF->prec_ints[tmpc*istep + k*jstep + tmpa*kstep + i]
                                          -postHF->prec_ints[   i*istep + k*jstep + tmpa*kstep + tmpc])*CC->T1[k*nvir+c];
            }

            for(int c=0;c<nvir;c++)
                for(int d=0;d<nvir;d++)
                   for(int k=0;k<nocc;k++){
                       tmpc = c+nocc;
                       tmpa = a+nocc;
                       tmpd = d+nocc;
                       //r[ia] += [2*(ckad)-(dkac)]tau[kicd]
                       r1[i*nvir +a] += (2*postHF->prec_ints[tmpc*istep + k*jstep + tmpa*kstep + tmpd]
                                          -postHF->prec_ints[tmpd*istep + k*jstep + tmpa*kstep + tmpc])*
                                                         tau[k*Ti_step + i*Tj_step + c*Ta_step + d];
                   }

            for(int k=0;k<nocc;k++)
                for(int l=0;l<nocc;l++)
                    for(int c=0;c<nvir;c++){
                        tmpc = c+nocc;
                        //r1[ia] += [2(ckil)-(clik)]*tau[klca]
                        r1[i*nvir +a] -= (2*postHF->prec_ints[tmpc*istep + k*jstep + i*kstep + l]
                                           -postHF->prec_ints[tmpc*istep + l*jstep + i*kstep + k])*
                                                          tau[k*Ti_step + l*Tj_step + c*Ta_step + a];

                    }
            r1[i*nvir +a]  /= ( f[i*nroao+i] -f[(nocc+a)*nroao +nocc+a]);
            }

    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int a=0;a<nvir;a++)
                for(int b=0;b<nvir;b++){

    for(int c=0;c<nvir;c++){
        //P[ijab] += Gca[ca] * T2[ijcb]
        Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] += Gca[c*nvir +a]*CC->T2[i*Ti_step + j*Tj_step + c*Ta_step + b];
    }

    for(int k=0;k<nocc;k++){
        //P[ijab] -= Gik[ik] * T2[ijab]
        Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -= Gik[i*nocc +k]*CC->T2[k*Ti_step + j*Tj_step + a*Ta_step + b];
    }

    for(int c=0;c<nvir;c++)
      for(int k=0;k<nocc;k++){
          tmpc = c + nocc;
          tmpb = b + nocc;
          tmpa = a + nocc;
          //P[ijab] += [(iabc)-(ikbc)*T1[ka]]T1[jc]
          Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] += (postHF->prec_ints[i*istep + tmpa*jstep + tmpb*kstep + tmpc]
                                                    -postHF->prec_ints[i*istep + k*jstep    + tmpb*kstep + tmpc]
                                                     *CC->T1[k*nvir+a])
                                                     *CC->T1[j*nvir+c];
      }
    for(int c=0;c<nvir;c++)
      for(int k=0;k<nocc;k++){
          //P[ijab] -= [(aijk)+(aick)*T1[jc]]*T1[kb]
          tmpc = c + nocc;
          Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -= (postHF->prec_ints[tmpa*istep + i*jstep + j*kstep    + k]
                                                    +postHF->prec_ints[tmpa*istep + i*jstep + tmpc*kstep + k]
                                                    *CC->T1[j*nvir+c])*CC->T1[k*nvir+b];
    }
    for(int c=0;c<nvir;c++)
      for(int k=0;k<nocc;k++){
          //P[ijab] +=(J[icak]-0.5K[icka])(2*T2[kjcb]-T2[kjbc])
          Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] +=      (Jicak[i*J_istep + c*J_cstep + a*J_astep + k]
                                                    - 0.5*Kicka[i*K_istep + c*K_cstep + k*K_kstep + a])*
                                                    (2*CC->T2[k*Ti_step + j*Tj_step + c*Ta_step + b]
                                                     - CC->T2[k*Ti_step + j*Tj_step + b*Ta_step + c]);
    }
    for(int c=0;c<nvir;c++)
      for(int k=0;k<nocc;k++){

          //P[ijab] -= 1/2(K[icka]*T2[kjbc]
          Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -= 0.5*(Kicka[i*K_istep + c*K_cstep + k*K_kstep + a]
                                                       *CC->T2[k*Ti_step + j*Tj_step + b*Ta_step + c]);
      }
    for(int c=0;c<nvir;c++)
      for(int k=0;k<nocc;k++){
          //P[ijab] -= (K[ickb]*T2[kjac]
          Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -=     (Kicka[i*K_istep + c*K_cstep + k*K_kstep + b]
                                                       *CC->T2[k*Ti_step + j*Tj_step + a*Ta_step + c]);
      }


 }

    for(int i=0;i<nocc;i++)
        for(int j=0;j<nocc;j++)
            for(int a=0;a<nvir;a++)
               for(int b=0;b<nvir;b++){
                   tmpa = a+nocc;
                   tmpb = b+nocc;
                   //r2[ijab] = (iajb)
                   r2[i*Ti_step+j*Tj_step+a*Ta_step+b] = postHF->prec_ints[i*istep + tmpa*jstep + j*kstep + tmpb];
                   for(int k=0;k<nocc;k++)
                       for(int l=0;l<nocc;l++){
                          //r2[ijab] = A[ijkl]*klab
                          r2[i*Ti_step+j*Tj_step+a*Ta_step+b] += Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l]
                                                                 * tau[k*Ti_step + l*Tj_step + a*Ta_step + b];
                       }

                   for(int c=0;c<nvir;c++)
                       for(int d=0;d<nvir;d++){
                            //r2[ijab] =B[cdab]*tau[ijcd]
                            r2[i*Ti_step+j*Tj_step+a*Ta_step+b] += Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b]
                                                                   * tau[i*Ti_step + j*Tj_step + c*Ta_step + d];
                       }


                   r2[i*Ti_step+j*Tj_step+a*Ta_step+b] += Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] + Pijab[j*Ti_step+i*Tj_step+b*Ta_step+a];

                   r2[i*Ti_step+j*Tj_step+a*Ta_step+b] /=  (f[i*nroao+i]+f[j*nroao + j]-f[(nocc+a)*nroao+(nocc+a)]-f[(nocc+b)*nroao + nocc+b]);


        }




    double* swapper = CC->T2;
    CC->T2  = r2;
    r2 = swapper;
    swapper = CC->T1;
    CC->T1  = r1;
    r1 = swapper;

    std::cout<<"T1:\n";
    for (int i=0;i<nocc;i++)
        for (int a=0;a<nvir;a++){
            if (abs(CC->T1[i*nvir+a]) > 1E-5) {
                std::cout<<i<<" "<<a<<" "<< CC->T1[i*nvir+a]<<std::endl;
            }

        }


    std::cout<<"T2:\n";
    for (int i=0;i<nocc;i++)
        for (int j=0;j<nocc;j++)
            for (int a=0;a<nvir;a++)
                for (int b=0;b<nvir;b++){
                    if (abs(CC->T2[i*Ti_step+j*Tj_step+a*Ta_step+b]) > 1E-5){
                        std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" "<<CC->T2[i*Ti_step+j*Tj_step+a*Ta_step+b]<<"\n";
                    }
                }


}







}







