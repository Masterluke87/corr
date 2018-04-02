#include "../ccsd.h"
#include "ccsd_intermediates_restr.h"
#include <cblas.h>
#include <lapacke.h>


void dgels_( char*, int*, int*, int*, double*, int*,
             double*, int*, double*, int*, int* );


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

void ccsd_guess_restr(cc_helper *CC,pHF* postHF){

    int nvir =  CC->nvir;
    int nocc =  CC->nocc;
    int nroao = CC->nroao;

    double* prec_ints = postHF->prec_ints;
    double* f = CC->f;

    double* T1 = CC->T1;
    double* T2 = CC->T2;



    memset(T1,0,nocc*nvir*sizeof(double));
    memset(T2,0,nocc*nocc*nvir*nvir*sizeof(double));

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
                    T2[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] = -(postHF->prec_ints[Ti*nistep + tmpA*njstep + Tj*nkstep + tmpB])/
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
        for(int a=0; a<nvir; a++)
        {
            tmpa=a+nocc;
            sum+= f[ i * nroao + tmpa ]*T1[ i * nvir + a ];
        }
    sum*= 2.0;
    CC->Ecc+=sum;

    //\frac{t^{dc}_{kl} v^{kl}_{dc}}{4}
    sum = 0.0;
    for(int a=0; a < nvir; a++)
        for(int b=0; b < nvir; b++)
            for(int i=0; i < nocc; i++)
                for(int j=0; j < nocc; j++)
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

void allocate_amplitudes_intermediates_restr(cc_helper* CC,cc_intermediates_restr* CC_int){
    int nocc=CC->nocc;
    int nvir=CC->nvir;

    long int mem_a = 2*(nocc*nvir) +                            //T1
                     2*(nocc*nocc*nvir*nvir);                   //T2

    long int mem_int = (nvir*nvir) +                            //Hca
                       (nvir*nvir) +                            //Gca
                       (nocc*nocc) +                            //Hik
                       (nocc*nocc) +                            //Gik
                       (nocc*nvir) +                            //Hck
                       (nocc*nocc*nocc*nocc)+                   //Aijkl
                       (nvir*nvir*nvir*nvir)+                   //Babcd
                       (nocc*nvir*nvir*nocc)+                   //Jicak
                       (nocc*nocc*nvir*nvir)+                   //Kicka
                       (nocc*nocc*nvir*nvir)+                    //Tau
                       (nocc*nocc*nvir*nvir);                   //Pijab

    std::cout<<"\nAllocating "<<mem_a/1024/1024 <<"MB("<<mem_a/1024/1024/1024<<"Gb) for amplitues\n";
    std::cout<<"Allocating "<<mem_int/1024/1024<<"MB("<<mem_int/1024/1024/1024<<"Gb) for intermediates\n";
    std::cout.flush();
    CC->pMem = new double[(mem_a+mem_int)];
    long int inc = 0;
    memset(CC->pMem,0,(mem_a+mem_int)*sizeof(double));

    CC->T1  = &(CC->pMem[inc]); inc+=nocc*nvir;
    CC->T1n = &(CC->pMem[inc]); inc+=nocc*nvir;
    CC->T2  = &(CC->pMem[inc]); inc+=nocc*nocc*nvir*nvir;
    CC->T2n = &(CC->pMem[inc]); inc+=nocc*nocc*nvir*nvir;

    CC_int->Hca = &(CC->pMem[inc]); inc+=nvir*nvir;
    CC_int->Gca = &(CC->pMem[inc]); inc+=nvir*nvir;

    CC_int->Gik = &(CC->pMem[inc]); inc+=nocc*nocc;
    CC_int->Hik = &(CC->pMem[inc]); inc+=nocc*nocc;

    CC_int->Hck = &(CC->pMem[inc]); inc+=nocc*nvir;

    CC_int->Aijkl = &(CC->pMem[inc]); inc+=nocc*nocc*nocc*nocc;
    CC_int->Bcdab = &(CC->pMem[inc]); inc+=nvir*nvir*nvir*nvir;

    CC_int->Jicak = &(CC->pMem[inc]); inc+=nocc*nvir*nvir*nocc;
    CC_int->Kicka = &(CC->pMem[inc]); inc+=nocc*nvir*nocc*nvir;

    CC_int->tau   = &(CC->pMem[inc]); inc+=nvir*nvir*nocc*nocc;
    CC_int->Pijab = &(CC->pMem[inc]);
}


void ccsd_build_intermediates_restr(cc_helper *CC,cc_intermediates_restr* CC_int,pHF* postHF){
    build_tau(CC,CC_int,postHF);
    build_Hik(CC,CC_int,postHF);
    build_Hca(CC,CC_int,postHF);
    build_Hck(CC,CC_int,postHF);
    build_Gik(CC,CC_int,postHF);
    build_Gca(CC,CC_int,postHF);
    build_Aikl(CC,CC_int,postHF);
    build_Bcdab(CC,CC_int,postHF);
    build_Jicak(CC,CC_int,postHF);
    build_Kicka(CC,CC_int,postHF);

}

void update_T1_restr(cc_helper *CC,cc_intermediates_restr* CC_int,pHF* postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;



    double* r1  = CC->T1n;
    double* tau = CC_int->tau;
    double* Hck = CC_int->Hck;
    double* Hca = CC_int->Hca;
    double* Hik = CC_int->Hik;
    double* f   = CC->f;

    int tmpa,tmpc,tmpd;

    memset(r1,0,nocc*nvir*sizeof(double));
    for(int i=0; i<nocc; i++)
        for(int a=0; a<nvir; a++) {
            for(int c=0; c<nvir; c++) {
                //r[ia] += Hca[ca]T1[ic]
                r1[i*nvir +a] += Hca[c*nvir+a]*CC->T1[i*nvir+c];
            }
            for(int k=0; k<nocc; k++) {
                //r[ia] -= Hik[ik]T1[ka]
                r1[i*nvir +a] -= Hik[i*nocc+k]*CC->T1[k*nvir+a];
            }
            for(int c=0; c<nvir; c++)
                for(int k=0; k<nocc; k++) {
                    //r[ia]+= Hck[ck]*(2*T2[kica] - T2[ikca] + T1[ic]*T1[ka])
                    r1[i*nvir +a] +=  Hck[c*nocc+k]*
                                     (2*CC->T2[k*Ti_step + i*Tj_step + c*Ta_step + a]
                                      - CC->T2[i*Ti_step + k*Tj_step + c*Ta_step + a]
                                      + CC->T1[i*nvir+c]*CC->T1[k*nvir+a]);
                }
            for(int c=0; c<nvir; c++)
                for(int k=0; k<nocc; k++) {
                    tmpc = c+nocc;
                    tmpa = a+nocc;
                    // r[ia]+= [2(ckai)-(ikac)]T1[kc]
                    r1[i*nvir +a] += (2*postHF->prec_ints[tmpc*istep + k*jstep + tmpa*kstep + i]
                                      -postHF->prec_ints[   i*istep + k*jstep + tmpa*kstep + tmpc])*CC->T1[k*nvir+c];
                }

            for(int c=0; c<nvir; c++)
                for(int d=0; d<nvir; d++)
                    for(int k=0; k<nocc; k++) {
                        tmpc = c+nocc;
                        tmpa = a+nocc;
                        tmpd = d+nocc;
                        //r[ia] += [2*(ckad)-(dkac)]tau[kicd]
                        r1[i*nvir +a] += (2*postHF->prec_ints[tmpc*istep + k*jstep + tmpa*kstep + tmpd]
                                          -postHF->prec_ints[tmpd*istep + k*jstep + tmpa*kstep + tmpc])*
                                         tau[k*Ti_step + i*Tj_step + c*Ta_step + d];
                    }

            for(int k=0; k<nocc; k++)
                for(int l=0; l<nocc; l++)
                    for(int c=0; c<nvir; c++) {
                        tmpc = c+nocc;
                        //r1[ia] += [2(ckil)-(clik)]*tau[klca]
                        r1[i*nvir +a] -= (2*postHF->prec_ints[tmpc*istep + k*jstep + i*kstep + l]
                                          -postHF->prec_ints[tmpc*istep + l*jstep + i*kstep + k])*
                                         tau[k*Ti_step + l*Tj_step + c*Ta_step + a];

                    }
            r1[i*nvir +a]  /= ( f[i*nroao+i] -f[(nocc+a)*nroao +nocc+a]);
        }


}

void update_T2_restr(cc_helper *CC,cc_intermediates_restr* CC_int,pHF* postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* r2 = CC->T2n;
    double* tau = CC_int->tau;
    double* Aijkl = CC_int->Aijkl;
    double* Bcdab = CC_int->Bcdab;
    double* Pijab = CC_int->Pijab;
    double* f = CC->f;

    int A_kstep = nocc;
    int A_jstep = A_kstep*nocc;
    int A_istep = A_jstep*nocc;

    int B_astep = nvir;
    int B_dstep = B_astep*nvir;
    int B_cstep = B_dstep*nvir;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    memset(r2,0,nocc*nocc*nvir*nvir*sizeof(double));
    memset(Pijab,0,nocc*nocc*nvir*nvir*sizeof(double));



    build_Pijab(CC,CC_int,postHF);
    int tmpa,tmpb;
    for(int i=0; i<nocc; i++)
        for(int j=0; j<nocc; j++)
            for(int a=0; a<nvir; a++)
                for(int b=0; b<nvir; b++) {
                    tmpa = a+nocc;
                    tmpb = b+nocc;
                    //r2[ijab] = (iajb)
                    r2[i*Ti_step+j*Tj_step+a*Ta_step+b] = postHF->prec_ints[i*istep + tmpa*jstep + j*kstep + tmpb];
                    for(int k=0; k<nocc; k++)
                        for(int l=0; l<nocc; l++) {
                            //r2[ijab] = A[ijkl]*klab
                            r2[i*Ti_step+j*Tj_step+a*Ta_step+b] += Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l]
                                                                   * tau[k*Ti_step + l*Tj_step + a*Ta_step + b];
                        }

                    for(int c=0; c<nvir; c++)
                        for(int d=0; d<nvir; d++) {
                            //r2[ijab] =B[cdab]*tau[ijcd]
                            r2[i*Ti_step+j*Tj_step+a*Ta_step+b] += Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b]
                                                                   * tau[i*Ti_step + j*Tj_step + c*Ta_step + d];
                        }


                    r2[i*Ti_step+j*Tj_step+a*Ta_step+b] += Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] + Pijab[j*Ti_step+i*Tj_step+b*Ta_step+a];

                    r2[i*Ti_step+j*Tj_step+a*Ta_step+b] /=  (f[i*nroao+i]+f[j*nroao + j]-f[(nocc+a)*nroao+(nocc+a)]-f[(nocc+b)*nroao + nocc+b]);


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

    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;



    int diis_size = 5;
    int d_size = diis_size+1;
    double* ampl_mem = new double[(diis_size+1)*((nocc*nvir)+(nocc*nocc*nvir*nvir))];
    double* error    = new double[(diis_size)*((nocc*nvir)+(nocc*nocc*nvir*nvir))];

    memset(ampl_mem,0,(1+diis_size)*((nocc*nvir)+(nocc*nocc*nvir*nvir))*sizeof(double));
    std::vector<double*> ampl(diis_size+1);
    std::vector<double*> err(diis_size);

    for(int i=0; i<(diis_size+1); i++) {
        ampl[i] = &(ampl_mem[i*((nocc*nvir)+(nocc*nocc*nvir*nvir))]);
    }

    for(int i=0; i<(diis_size); i++) {
        err[i] = &(error[i*((nocc*nvir)+(nocc*nocc*nvir*nvir))]);
    }


    double* tmp0;
    bool use_diis = true;
    double* B  = new double[(diis_size+1)*(diis_size+1)];
    double* B2 = new double[(diis_size+1)*(diis_size+1)];
    double* c  = new double[(diis_size+1)];

    double maxel = 0.0;


    form_fock_restr(CC,postHF);
    allocate_amplitudes_intermediates_restr(CC,CC_int);
    calc_e_check_restr(CC,postHF);
    ccsd_guess_restr(CC,postHF);

    for(int iter=0; iter<20; iter++) {

        ccsd_energy_restr(CC,postHF);
        ccsd_build_intermediates_restr(CC,CC_int,postHF);

        update_T1_restr(CC,CC_int,postHF);
        update_T2_restr(CC,CC_int,postHF);

        //first move all up

        if (use_diis) {
            tmp0 = ampl[0];
            for(int i=1; i<(diis_size+1); i++)
                ampl[i-1]  = ampl[i];

            ampl.back() = tmp0;

            tmp0 = err[0];
            for(int i=1; i<(diis_size); i++)
                err[i-1]  = err[i];

            err.back() = tmp0;


            //add error vector



            for (int i=0; i<(nocc*nvir); i++)
                ampl.back()[i] = CC->T1n[i];
            for (int i=0; i<(nocc*nocc*nvir*nvir); i++)
                ampl.back()[i+nocc*nvir] = CC->T2n[i];

            for (int i=0; i<(nocc*nvir); i++)
                err.back()[i] = CC->T1n[i]-CC->T1[i];

//            for (int i=0; i<(nocc*nocc*nvir*nvir); i++)
//                err.back()[i+nocc*nvir] = CC->T2n[i]-CC->T2[i];


            for (int i=0; i<(nocc*nvir+nocc*nocc*nvir*nvir); i++)
               err.back()[i] = ampl.back()[i] - ampl[diis_size-1][i];


            for(int d=0; d<(diis_size+1); d++) {
                std::cout<<d<<" AMP:";
                for (int i=0; i<(nocc*nvir+nocc*nocc*nvir*nvir); i++) {
                    std::cout<<" "<<ampl[d][i];
                }
                std::cout<<"\n";
            }


            for(int d=0; d<diis_size; d++) {
                std::cout<<d<<" ERR:";
                for (int i=0; i<(nocc*nvir+nocc*nocc*nvir*nvir); i++) {
                    std::cout<<" "<<err[d][i];
                }
                std::cout<<"\n";
            }



            if (iter>=diis_size) {
                maxel =0.0;
                memset(B,0,d_size*d_size*sizeof(double));
                memset(B2,0,d_size*d_size*sizeof(double));
                for (int i=0; i<diis_size; i++) {
                    for (int j=0; j<diis_size; j++) {
                        B[i*d_size + j] += cblas_ddot((nocc*nvir)+(nocc*nocc*nvir*nvir),err[i],1,err[j],1);
                        for(int k=0; k<((nocc*nvir)+(nocc*nocc*nvir*nvir)); k++) {
                            B2[i*d_size + j] += (ampl[i+1][k]-ampl[i][k]) * (ampl[j+1][k]-ampl[j][k]);
                        }

                        if (maxel<fabs(B[i*d_size + j]))
                            maxel = B[i*d_size + j];
                    }
                    B[i*d_size + d_size-1] = -1;
                    B[(d_size-1)*d_size+i] = -1;
                    c[i] = 0.0;

                }
                //std::cout<<maxel<<"\n";
                B[d_size*d_size-1] = 0.0;
                c[d_size-1] = -1;
                for (int i=0; i<diis_size; i++) {
                    for (int j=0; j<diis_size; j++) {
                        B[i*d_size + j] /= maxel;
                        B2[i*d_size + j] /= maxel;
                    }
                }
/*
                for (int i=0; i<d_size; i++) {
                    for (int j=0; j<d_size; j++) {
                        std::cout<<B[i*d_size + j]<<" ";
                    }
                    std::cout<<std::endl;
                }
                std::cout.flush();
                for (int i=0; i<d_size; i++) {
                    for (int j=0; j<d_size; j++) {
                        std::cout<<B2[i*d_size + j]<<" ";
                    }
                    std::cout<<std::endl;
                }
                std::cout.flush();
*/
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

                memset(CC->T1n,0,nocc*nvir*sizeof(double));
                memset(CC->T2n,0,nocc*nocc*nvir*nvir*sizeof(double));

                for(int j=0; j<diis_size; j++) {
                    for(int i=0; i<(nocc*nvir); i++) {
                        CC->T1n[i] += c[j] * ampl[j+1][i];
                    }
                }
                for(int j=0; j<diis_size; j++) {
                    for(int i=0; i<(nocc*nocc*nvir*nvir); i++) {
                        CC->T2n[i] += c[j] * ampl[j+1][i+nocc*nvir];
                    }
                }
                for(int i=0;i<((nocc*nvir);i++){
                    std::cout<<CC->T1n[i];
                }
                    for(int i=0;i<((nocc*nvir);i++){
                        std::cout<<CC->T1n[i];
                    }



            }


        }

        double* swapper = CC->T2;
        CC->T2  = CC->T2n;
        CC->T2n = swapper;
        swapper = CC->T1;
        CC->T1  = CC->T1n;
        CC->T1n = swapper;

        /*
           std::cout<<"T1:\n";
           for (int i=0; i<nocc; i++)
            for (int a=0; a<nvir; a++) {
                if (abs(CC->T1[i*nvir+a]) > 1E-5) {
                    std::cout<<i<<" "<<a<<" "<< CC->T1[i*nvir+a]<<std::endl;
                }

            }


           std::cout<<"T2:\n";
           for (int i=0; i<nocc; i++)
            for (int j=0; j<nocc; j++)
                for (int a=0; a<nvir; a++)
                    for (int b=0; b<nvir; b++) {
                        if (abs(CC->T2[i*Ti_step+j*Tj_step+a*Ta_step+b]) > 1E-5) {
                            std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" "<<CC->T2[i*Ti_step+j*Tj_step+a*Ta_step+b]<<"\n";
                        }
                    }
         */

    }
}
