#include "ccsd_intermediates.h"
#include "ccsd.h"



void allocate_intermediates_memory(cc_helper *CC,cc_intermediates* CC_int){
    int nocc = CC->nocc;
    int nvir = CC->nvir;

    long int mem = (nvir*nvir) + (nocc*nocc) + (nocc*nvir) +
            (nocc*nocc*nocc*nocc)+
            (nvir*nvir*nvir*nvir)+
            (nocc*nvir*nvir*nocc)+
            (nocc*nocc*nvir*nvir)+
            (nocc*nocc*nvir*nvir);
    CC_int->pMem = new double[mem];
    long int inc = 0;

    CC_int->Fae = &(CC_int->pMem[inc]); inc+=nvir*nvir;
    CC_int->Fmi = &(CC_int->pMem[inc]); inc+=nocc*nocc;
    CC_int->Fme = &(CC_int->pMem[inc]); inc+=nocc*nvir;

    CC_int->Wmnij = &(CC_int->pMem[inc]); inc+=nocc*nocc*nocc*nocc;
    CC_int->Wabef = &(CC_int->pMem[inc]); inc+=nvir*nvir*nvir*nvir;
    CC_int->Wmbej = &(CC_int->pMem[inc]); inc+=nocc*nvir*nvir*nocc;

    CC_int->tau   = &(CC_int->pMem[inc]); inc+=nocc*nocc*nvir*nvir;
    CC_int->tau_s = &(CC_int->pMem[inc]);








}


void build_Fae(cc_helper* CC, cc_intermediates *CC_int, pHF* postHF){
    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int nroao = CC->nroao;
    int norb  = CC->norb;


    double* Fae   = CC_int->Fae;
    double* tau   = CC_int->tau;
    double* tau_s = CC_int->tau_s;
    double* prec_ints = postHF->prec_ints;
    double* T1 = CC->T1;
    double* f = CC->f;

    int tmpa,tmpe,tmpf;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

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

}

void build_Fmi(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF){
    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int norb = CC->norb;
    int nroao = CC->nroao;

    double* Fmi   = CC_int->Fmi;
    double* tau_s = CC_int->tau_s;
    double* f     = CC->f;
    double* T1    = CC->T1;

    double* prec_ints = postHF->prec_ints;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int tmpe,tmpf;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;


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
}
void build_Fme(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF){
    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int norb = CC->norb;
    int nroao = CC->nroao;

    double* Fme = CC_int->Fme;
    double* f   = CC->f;
    double* T1  = CC->T1;
    double* prec_ints = postHF->prec_ints;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;



    int tmpe,tmpf;

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

}
void build_Wmnij(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF){

    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int norb = CC->norb;
    int nroao = CC->nroao;

    double* Wmnij = CC_int->Wmnij;
    double* prec_ints = postHF->prec_ints;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;


    double* T1 = CC->T1;
    double* tau = CC_int->tau;

    int tmpe,tmpf;

    int m_step = nocc*nocc*nocc;
    int n_step = nocc*nocc;
    int i_step = nocc;


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

}
void build_Wabef(cc_helper* CC,cc_intermediates* CC_int,pHF* postHF){

    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int norb = CC->norb;
    int nroao = CC->nroao;

    double* Wabef = CC_int->Wabef;
    double* T1 = CC->T1;
    double* tau = CC_int->tau;

    double* prec_ints = postHF->prec_ints;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int a_step = nvir*nvir*nvir;
    int b_step = nvir*nvir;
    int e_step = nvir;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;
    int tmpe,tmpf,tmpa,tmpb;
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

}

void build_Wmbej(cc_helper* CC, cc_intermediates *CC_int, pHF* postHF){
    int nvir = CC->nvir;
    int nocc = CC->nocc;
    int norb = CC->norb;
    int nroao = CC->nroao;

    double* Wmbej = CC_int->Wmbej;
    double* T1 = CC->T1;
    double* T2 = CC->T2;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    double* prec_ints = postHF->prec_ints;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int m_step = nvir*nvir*nocc;
    int b_step = nvir*nocc;
    int e_step = nocc;

    int tmpe,tmpf,tmpa,tmpb;


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

}
