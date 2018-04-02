#include "ccsd_intermediates_restr.h"



void build_tau(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir = CC->nvir;
    int nocc = CC->nocc;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;
    for(int i=0; i<nocc; i++)
        for(int j=0; j<nocc; j++)
            for(int a=0; a<nvir; a++)
                for(int b=0; b<nvir; b++) {
                    //tau[ijab] = t2[ijab]+t1[ia]*t1[jb]
                    CC_int->tau[i*Ti_step+j*Tj_step+a*Ta_step+b] = CC->T2[i*Ti_step+j*Tj_step+a*Ta_step+b] + CC->T1[i*nvir+a]* CC->T1[j*nvir+b];
                }

}

void build_Hik(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* Hik = CC_int->Hik;
    double* tau = CC_int->tau;


    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    int tmpc,tmpd;
    for(int i=0; i<nocc; i++)
        for(int k=0; k<nocc; k++) {
            Hik[i*nocc+k] = 0.0;
            for(int l=0; l<nocc; l++)
                for(int c=0; c<nvir; c++)
                    for(int d=0; d<nvir; d++) {
                        tmpc = c+nocc;
                        tmpd = d+nocc;
                        //Hik[ik] = [2(kcld)-(kdlc)]tau[ilcd]
                        Hik[i*nocc+k] += (2*postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd] - postHF->prec_ints[k*istep + tmpd*jstep + l*kstep + tmpc])
                                         *tau[i*Ti_step + l*Tj_step +  c*Ta_step +  d];
                    }

        }

}

void build_Gik(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* Hik = CC_int->Hik;
    double* Gik = CC_int->Gik;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    int tmpc = 0;
    for (int i=0; i<nocc; i++)
        for(int k=0; k<nocc; k++) {
            Gik[i*nocc + k] = Hik[i*nocc+k];
            for(int c=0; c<nvir; c++)
                for(int l=0; l<nocc; l++) {
                    tmpc = c+nocc;
                    //Gik[ik]=[2*(ikcl) - (ilck)]*T1[lc]
                    Gik[i*nocc + k] += (2*postHF->prec_ints[i*istep + k*jstep + tmpc*kstep + l]
                                        - postHF->prec_ints[i*istep + l*jstep + tmpc*kstep + k])*CC->T1[l*nvir+c];
                }
        }

}

void build_Hca(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* Hca = CC_int->Hca;
    double* tau = CC_int->tau;

    double sum = 0.0;
    int tmpc = 0;
    int tmpd;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    for(int a=0; a<nvir; a++)
        for(int c=0; c<nvir; c++) {
            Hca[c*nvir+a] = 0.0;
            sum = 0.0;

            for(int k=0; k<nocc; k++)
                for(int l=0; l<nocc; l++)
                    for(int d=0; d<nvir; d++) {
                        tmpc = c+nocc;
                        tmpd = d+nocc;
                        //[2*(kcld) - (kdlc)]tau[klad]
                        sum += (2*postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd]
                                - postHF->prec_ints[k*istep + tmpd*jstep + l*kstep + tmpc])
                               *tau[k*Ti_step + l*Tj_step +  a*Ta_step +  d];

                    }

            Hca[c*nvir+a] =-sum;
        }
}

void build_Gca(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    double* Gca = CC_int->Gca;
    double* Hca = CC_int->Hca;


    int tmpa,tmpc,tmpd;
    for (int a=0; a<nvir; a++)
        for(int c=0; c<nvir; c++) {
            Gca[c*nvir+a] = Hca[c*nvir+a];
            for(int d=0; d<nvir; d++)
                for(int k=0; k<nocc; k++) {
                    tmpa = a+nocc;
                    tmpc = c+nocc;
                    tmpd = d+nocc;
                    //Gca[ca] += [2(acdk)- (adck)]T1[kd]
                    Gca[c*nvir+a] += (2*postHF->prec_ints[tmpa*istep + tmpc*jstep + tmpd*kstep + k]
                                      - postHF->prec_ints[tmpa*istep + tmpd*jstep + tmpc*kstep + k])*CC->T1[k*nvir+d];
                }
        }


}

void build_Hck(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* Hck = CC_int->Hck;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;


    int tmpc;
    int tmpd;

    for(int c=0; c<nvir; c++)
        for(int k=0; k<nocc; k++) {
            Hck[c*nocc+k] = 0.0;
            for(int l=0; l<nocc; l++)
                for(int d=0; d<nvir; d++) {
                    tmpc = c+nocc;
                    tmpd = d+nocc;
                    //H[ck] = [2*(kcld)-(kdlc)]T1[ld]
                    Hck[c*nocc+k]+= (2*postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd]
                                     - postHF->prec_ints[k*istep + tmpd*jstep + l*kstep + tmpc])*CC->T1[l*nvir+d];

                }
        }
}

void build_Aikl(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* Aijkl = CC_int->Aijkl;
    double* tau   = CC_int->tau;
    int A_kstep = nocc;
    int A_jstep = A_kstep*nocc;
    int A_istep = A_jstep*nocc;
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    int tmpc;
    int tmpd;
    for(int i=0; i<nocc; i++)
        for(int j=0; j<nocc; j++)
            for(int k=0; k<nocc; k++)
                for(int l=0; l<nocc; l++) {
                    //A[ijkl] = (ikjl)
                    Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] = postHF->prec_ints[i*istep + k*jstep + j*kstep + l];
                    for(int c=0; c<nvir; c++) {
                        tmpc = c+nocc;
                        //A[ijkl] += (ikcl)T[jc]
                        Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] += postHF->prec_ints[i*istep + k*jstep + tmpc*kstep + l]*CC->T1[j*nvir+c];
                    }
                    for(int c=0; c<nvir; c++) {
                        tmpc = c + nocc;
                        //A[ijkl] += (ckjl)T[ic]
                        Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] += postHF->prec_ints[tmpc*istep + k*jstep + j*kstep + l]*CC->T1[i*nvir+c];
                    }
                    for(int c=0; c<nvir; c++)
                        for(int d=0; d<nvir; d++) {
                            tmpc = c + nocc;
                            tmpd = d + nocc;
                            //A[ijkl]+= (kcld)*tau[ijcd]
                            Aijkl[i*A_istep + j*A_jstep + k*A_kstep + l] += postHF->prec_ints[k*istep + tmpc*jstep + l*kstep + tmpd]*tau[i*Ti_step + j*Tj_step + c*Ta_step +d];

                        }
                }
}

void build_Bcdab(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* Bcdab = CC_int->Bcdab;

    int B_astep = nvir;
    int B_dstep = B_astep*nvir;
    int B_cstep = B_dstep*nvir;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int tmpa,tmpb,tmpc,tmpd;
    for(int c=0; c<nvir; c++)
        for(int d=0; d<nvir; d++)
            for(int a=0; a<nvir; a++)
                for(int b=0; b<nvir; b++) {
                    tmpa = a + nocc;
                    tmpb = b + nocc;
                    tmpc = c + nocc;
                    tmpd = d + nocc;
                    //B[cdab] = (acbd)
                    Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b] = postHF->prec_ints[tmpa*istep + tmpc*jstep + tmpb*kstep + tmpd];
                    for (int k=0; k<nocc; k++) {
                        //B[cdab] -= (acdk)T[kb]
                        Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b] -=  postHF->prec_ints[tmpa*istep + tmpc*jstep + tmpd*kstep + k]*CC->T1[k*nvir +b];  //(acdk)t1[kb]
                        //B[cdab] -= (bdck)T[ka]
                        Bcdab[c*B_cstep + d*B_dstep + a*B_astep +b] -=  postHF->prec_ints[tmpb*istep + tmpd*jstep + tmpc*kstep + k]*CC->T1[k*nvir +a];  //(bdck)t1[ka]
                    }



                }

}

void build_Jicak(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;

    double* Jicak = CC_int->Jicak;

    int tmpa,tmpb,tmpc,tmpd;

    int J_astep = nocc;
    int J_cstep = J_astep*nvir;
    int J_istep = J_cstep*nvir;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;



    for(int i=0; i<nocc; i++)
        for(int c=0; c<nvir; c++)
            for(int a=0; a<nvir; a++)
                for(int k=0; k<nocc; k++) {
                    tmpc = c+nocc;
                    tmpa = a+nocc;
                    //J[icak] = (aick)
                    Jicak[i*J_istep + c*J_cstep + a*J_astep +k] = postHF->prec_ints[tmpa*istep + i*jstep + tmpc*kstep + k];

                    for (int l=0; l<nocc; l++) {
                        //J[icak] -= (ilck)T[la]
                        Jicak[i*J_istep + c*J_cstep + a*J_astep +k] -= postHF->prec_ints[i*istep + l*jstep + tmpc*kstep + k]*
                                                                       CC->T1[l*nvir+a];
                    }

                    for (int d=0; d<nvir; d++) {
                        tmpd = d + nocc;
                        //J[icak] += (adck)T[id}
                        Jicak[i*J_istep + c*J_cstep + a*J_astep +k] += postHF->prec_ints[tmpa*istep + tmpd*jstep + tmpc*kstep + k]
                                                                       *CC->T1[i*nvir+d];
                    }

                    for (int d=0; d<nvir; d++)
                        for (int l=0; l<nocc; l++) {
                            tmpd = d + nocc;
                            //J[icak] -= 1/2(ckdl)(t2[ilda]+2*t[id]*t[la]
                            Jicak[i*J_istep + c*J_cstep + a*J_astep +k] -= 0.5*postHF->prec_ints[tmpc*istep + k*jstep + tmpd*kstep + l]*
                                                                           (CC->T2[i*Ti_step + l*Tj_step + d*Ta_step +a] +
                                                                            2*CC->T1[i*nvir+d] * CC->T1[l*nvir+a]);
                        }

                    for (int d=0; d<nvir; d++)
                        for (int l=0; l<nocc; l++) {
                            tmpd = d + nocc;
                            //J[icak] += 1/2[2(ckdl)-(dkcl)]t2[ilad]
                            Jicak[i*J_istep + c*J_cstep + a*J_astep +k] += 0.5*(2*postHF->prec_ints[tmpc*istep + k*jstep + tmpd*kstep + l]-
                                                                                postHF->prec_ints[tmpd*istep + k*jstep + tmpc*kstep + l])*CC->T2[i*Ti_step + l*Tj_step + a*Ta_step +d];


                        }


                }


}

void build_Kicka(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{
    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;


    int K_kstep = nvir;
    int K_cstep = K_kstep*nocc;
    int K_istep = K_cstep*nvir;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    double* Kicka = CC_int->Kicka;
    int tmpa,tmpc,tmpd;
    for(int i=0; i<nocc; i++)
        for(int c=0; c<nvir; c++)
            for(int k=0; k<nocc; k++)
                for(int a=0; a<nvir; a++)
                {
                    tmpa = a + nocc;
                    tmpc = c + nocc;
                    //K[icka] = (ikac)
                    Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] = postHF->prec_ints[i*istep + k*jstep + tmpa*kstep + tmpc];
                    for(int l=0; l<nocc; l++) {
                        //K[icka] -= (ikcl)T1[la]
                        Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] -= postHF->prec_ints[i*istep + k*jstep + tmpc*kstep + l]*CC->T1[l*nvir+a];
                    }
                    for(int d=0; d<nvir; d++) {
                        // K[icka]+= (dkac)*T1[id]
                        tmpd = d+nocc;
                        Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] += postHF->prec_ints[tmpd*istep + k*jstep + tmpa*kstep + tmpc]*CC->T1[i*nvir+d];
                    }
                    for (int d=0; d<nvir; d++)
                        for (int l=0; l<nocc; l++) {
                            tmpd = d + nocc;
                            //K[icka] -= 1/2(dkcl)(t2[ilda]*2T1[id]T1[la]
                            Kicka[i*K_istep + c*K_cstep + k*K_kstep + a] -= 0.5*postHF->prec_ints[tmpd*istep + k*jstep + tmpc*kstep + l]*
                                                                            (CC->T2[i*Ti_step + l*Tj_step + d*Ta_step+ a]
                                                                             + 2*CC->T1[i*nvir+d]*CC->T1[l*nvir+a]);


                        }

                }

}

void build_Pijab(cc_helper *CC, cc_intermediates_restr *CC_int, pHF *postHF)
{

    int nvir  = CC->nvir;
    int nocc  = CC->nocc;
    int nroao = CC->nroao;



    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    int J_astep = nocc;
    int J_cstep = J_astep*nvir;
    int J_istep = J_cstep*nvir;

    int K_kstep = nvir;
    int K_cstep = K_kstep*nocc;
    int K_istep = K_cstep*nvir;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    double* Pijab = CC_int->Pijab;
    double* Kicka = CC_int->Kicka;
    double* Jicak = CC_int->Jicak;
    double* Gca   = CC_int->Gca;
    double* Gik   = CC_int->Gik;





    int tmpa,tmpc,tmpb;

    for(int i=0; i<nocc; i++)
        for(int j=0; j<nocc; j++)
            for(int a=0; a<nvir; a++)
                for(int b=0; b<nvir; b++) {

                    for(int c=0; c<nvir; c++) {
                        //P[ijab] += Gca[ca] * T2[ijcb]
                        Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] += Gca[c*nvir +a]*CC->T2[i*Ti_step + j*Tj_step + c*Ta_step + b];
                    }

                    for(int k=0; k<nocc; k++) {
                        //P[ijab] -= Gik[ik] * T2[ijab]
                        Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -= Gik[i*nocc +k]*CC->T2[k*Ti_step + j*Tj_step + a*Ta_step + b];
                    }

                    for(int c=0; c<nvir; c++) {
                        tmpc = c + nocc;
                        tmpb = b + nocc;
                        tmpa = a + nocc;
                        Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] += postHF->prec_ints[i*istep + tmpa*jstep + tmpb*kstep + tmpc] *CC->T1[j*nvir+c];
                        for(int k=0; k<nocc; k++) {
                            Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] += -postHF->prec_ints[i*istep + k*jstep    + tmpb*kstep + tmpc]*CC->T1[k*nvir+a]*CC->T1[j*nvir+c];
                        }
                    }
                    double s1 = 0.0;
                    for(int k=0; k<nocc; k++) {
                        s1 = 0.0;
                        for(int c=0; c<nvir; c++) {
                            tmpc = c+ nocc;
                            s1 += postHF->prec_ints[tmpa*istep + i*jstep + tmpc*kstep + k]*CC->T1[j*nvir+c];
                        }
                        //P[ijab] -= [(aijk)+(aick)*T1[jc]]*T1[kb]
                        Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -= (postHF->prec_ints[tmpa*istep + i*jstep + j*kstep    + k] + s1
                                                                   )*CC->T1[k*nvir+b];
                    }

                    for(int c=0; c<nvir; c++)
                        for(int k=0; k<nocc; k++) {
                            //P[ijab] +=(J[icak]-0.5K[icka])(2*T2[kjcb]-T2[kjbc])
                            Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] +=      (Jicak[i*J_istep + c*J_cstep + a*J_astep + k]
                                                                            - 0.5*Kicka[i*K_istep + c*K_cstep + k*K_kstep + a])*
                                                                      (2*CC->T2[k*Ti_step + j*Tj_step + c*Ta_step + b]
                                                                       - CC->T2[k*Ti_step + j*Tj_step + b*Ta_step + c]);
                        }
                    for(int c=0; c<nvir; c++)
                        for(int k=0; k<nocc; k++) {

                            //P[ijab] -= 1/2(K[icka]*T2[kjbc]
                            Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -= 0.5*(Kicka[i*K_istep + c*K_cstep + k*K_kstep + a]
                                                                           *CC->T2[k*Ti_step + j*Tj_step + b*Ta_step + c]);
                        }
                    for(int c=0; c<nvir; c++)
                        for(int k=0; k<nocc; k++) {
                            //P[ijab] -= (K[ickb]*T2[kjac]
                            Pijab[i*Ti_step+j*Tj_step+a*Ta_step+b] -=     (Kicka[i*K_istep + c*K_cstep + k*K_kstep + b]
                                                                           *CC->T2[k*Ti_step + j*Tj_step + a*Ta_step + c]);
                        }


                }

}
