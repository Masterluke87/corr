#include "ccsd.h"
#include "../UTIL/util.h"
#include "ccsd_intermediates.h"
#include <cblas.h>

double get_integral(double* ints,long long int &istep,long long int &jstep,long long int &kstep,int  &i,int  &j,int &k,int &l){
    return ints[(i/2)*istep + (k/2)*jstep + (j/2)*kstep + l/2]*(i%2==k%2)*(j%2==l%2) -
           ints[(i/2)*istep + (l/2)*jstep + (j/2)*kstep + k/2]*(i%2==l%2)*(j%2==k%2);
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

    std::cout << "\nCCSD:\n----" << '\n';


    calc_ints_so(CC,postHF);
    allocate_amplitudes_intermediates(CC,CC_int);
    form_fock(CC,postHF);
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







    //Intermediates:


    /*
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



    //f^{k}_{c} t^{c}_{k}

    //Tau intermediates

    //Fae

    //Fmi


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



    }
*/
}
