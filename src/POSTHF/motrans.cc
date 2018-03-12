#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "motrans.h"


void MOtrans_mod_sym(systeminfo *sysinfo,OEints* onemats,TEints *twomats,pHF* postHF){
    double* MOs              = onemats->MOs;
    int nroao                = sysinfo->nroao;
    long long int nrofint    = sysinfo->nrofint;
    long long int* sortcount = twomats->sortcount;
    double* intval           = twomats->intval;
    unsigned short* intnums  = twomats->intnums;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;


}


void MOtrans_mod(systeminfo *sysinfo,OEints* onemats,TEints *twomats,pHF* postHF){
    double* MOs              = onemats->MOs;
    int nroao                = sysinfo->nroao;
    long long int nrofint    = sysinfo->nrofint;
    long long int* sortcount = twomats->sortcount;
    double* intval           = twomats->intval;
    unsigned short* intnums  = twomats->intnums;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

    /* one index at a time algo */
    long long int prec_mem = (long long int) nroao*(long long int) nroao* (long long int) nroao* (long long int) nroao;
    // (ij|kl) <- (ab|cd)
    double* I_ibcd   = new double[nroao*nroao*nroao*nroao];
    int off_i;

    std::cout << "\nMOTRANS:\n-------" << '\n';
    std::cout << "MOtrans: Need "<<2*prec_mem/1024/1024<<" MB\n";
    std::cout << "AO-Integrals will be deleted, frees " <<sysinfo->nrofint/1024/1024<<"MB\n";

    for(int i=0;i<nroao;i++){
        off_i = i*nroao;
        for(long long int x = 0; x < sortcount[0]; x++) {
            //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
             I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+2]*kstep+intnums[x*4+3]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)

        }
        //Cai(aa|cc) & Cci(cc|aa)
        for(long long int x = sortcount[0]; x < sortcount[1]; x++) {
            //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
             I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+2]*kstep+intnums[x*4+3]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
             I_ibcd[i*istep+intnums[x*4+2]*jstep+intnums[x*4+0]*kstep+intnums[x*4+1]] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
        }
        //  (ab|ab) |  C-0 | 123
        //  (ab|ba) |  C-0 | 132
        //  (ba|ab) |  C-1 | 023
        //  (ba|ba) |  C-1 | 032
        for(long long int x = sortcount[1]; x < sortcount[2]; x++) {
            //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
            I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+2]*kstep+intnums[x*4+3]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
            I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+3]*kstep+intnums[x*4+2]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc)
            I_ibcd[i*istep+intnums[x*4+0]*jstep+intnums[x*4+2]*kstep+intnums[x*4+3]] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(3)    (ba|cd)
            I_ibcd[i*istep+intnums[x*4+0]*jstep+intnums[x*4+3]*kstep+intnums[x*4+2]] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(4)    (ba|dc)
        }
        // (aa|cd) | C-0 | 0123
        // (aa|dc) | C-0 | 0132
        // (cd|aa) | C-2 | 2301
        // (dc|aa) | C-3 | 3201
        for(long long int x = sortcount[2]; x < sortcount[3]; x++) {
            //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
            I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+2]*kstep+intnums[x*4+3]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
            I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+3]*kstep+intnums[x*4+2]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc)
            I_ibcd[i*istep+intnums[x*4+3]*jstep+intnums[x*4+0]*kstep+intnums[x*4+1]] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
            I_ibcd[i*istep+intnums[x*4+2]*jstep+intnums[x*4+0]*kstep+intnums[x*4+1]] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(6)    (dc|ab)
        }
        // (ab|cd) | C-0 | 0123
        // (ab|dc) | C-0 | 0132
        // (ba|cd) | C-1 | 1023
        // (ba|dc) | C-1 | 1032
        // (cd|ab) | C-2 | 2301
        // (cd|ba) | C-2 | 2310
        // (dc|ab) | C-3 | 3201
        // (dc|ba) | C-3 | 3210

        for(long long int x = sortcount[3]; x < nrofint; x++) {
            //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
            I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+2]*kstep+intnums[x*4+3]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
            I_ibcd[i*istep+intnums[x*4+1]*jstep+intnums[x*4+3]*kstep+intnums[x*4+2]] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc)
            I_ibcd[i*istep+intnums[x*4+0]*jstep+intnums[x*4+2]*kstep+intnums[x*4+3]] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(3)    (ba|cd)
            I_ibcd[i*istep+intnums[x*4+0]*jstep+intnums[x*4+3]*kstep+intnums[x*4+2]] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(4)    (ba|dc)
            I_ibcd[i*istep+intnums[x*4+3]*jstep+intnums[x*4+0]*kstep+intnums[x*4+1]] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
            I_ibcd[i*istep+intnums[x*4+3]*jstep+intnums[x*4+1]*kstep+intnums[x*4+0]] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(6)    (dc|ab)
            I_ibcd[i*istep+intnums[x*4+2]*jstep+intnums[x*4+0]*kstep+intnums[x*4+1]] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(7)    (cd|ba)
            I_ibcd[i*istep+intnums[x*4+2]*jstep+intnums[x*4+1]*kstep+intnums[x*4+0]] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(8)    (dc|ba)
        }

    }
    std::cout << "Done 1st index..\n";


    //remove the two electron integrals
    delete[] twomats->intval;
    delete[] twomats->intnums;
    delete   twomats;
    //Allocate final space

    postHF->prec_ints =  new double[prec_mem];
    std::memset(postHF->prec_ints,0,prec_mem*sizeof(double));
    //2nd index
    //(ij|cd) <- C_b*(i|bcd)
    for(int i=0;i<nroao;i++){
        for(int j=0;j<nroao;j++){
          for(int c=0;c<nroao;c++){
              for(int d=0;d<nroao;d++){
                    for(int b=0;b<nroao;b++)
                        {
                        postHF->prec_ints[i*istep + j*jstep + c*kstep + d] += MOs[j*nroao+b]*I_ibcd[i*istep + b*jstep + c*kstep + d];
                        }
                }
            }
        }
    }
    std::cout << "Done 2nd index..\n";

    double* I_ijkd = I_ibcd;
    std::memset(I_ijkd,0,prec_mem*sizeof(double));
    //3rd index
    //(ij|kd) <- C_c(ij|cd)
    for(int i=0;i<nroao;i++){
        for(int j=0;j<nroao;j++){
            for(int k=0;k<nroao;k++){
                for(int d=0;d<nroao;d++){
                    for(int c=0;c<nroao;c++){
                        I_ijkd[i*istep + j*jstep + k*kstep + d] += MOs[k*nroao+c]*postHF->prec_ints[i*istep + j*jstep + c*kstep + d];
                    }
                }
            }
        }
    }
    std::cout << "Done 3rd index..\n";
    std::memset(postHF->prec_ints,0,nroao*nroao*nroao*nroao*sizeof(double));
    //(ij|kl) <- C_d(ij|kd)
    for(int i=0;i<nroao;i++){
        for(int j=0;j<nroao;j++){
            for(int k=0;k<nroao;k++){
                for(int l=0;l<nroao;l++){
                    for(int d=0;d<nroao;d++){
                       postHF->prec_ints[i*istep + j*jstep + k*kstep + l] += MOs[l*nroao+d]*I_ijkd[i*istep + j*jstep + k*kstep + d];
                    }
                }
            }
        }
    }
    std::cout << "Done 4th index..\n";

    delete[] I_ijkd;

}






void build_FMo(systeminfo *sysinfo,OEints* onemats, pHF* postHF)
    {


    int nroao    = sysinfo->nroao;
    double* Fmat = onemats->Fmat;

    postHF->FMo = new double[nroao*nroao];


    double* FMo  = postHF->FMo;
    double* MOs  = onemats->MOs;

    memset(FMo,0,nroao*nroao*sizeof(double));

    for (int i = 0; i < nroao; i++) {
        for (int j = 0; j < nroao; j++)
            for (int mu = 0; mu < nroao; mu++) {
                for (int nu = 0; nu < nroao; nu++) {
                    FMo[i*nroao + j] += MOs[i*nroao + mu] * Fmat[mu*nroao + nu] * MOs[j*nroao + nu];

            }
        }
    }
}


void MOtrans(systeminfo *sysinfo,OEints* onemats,TEints *twomats,pHF* postHF)
{
    double* MOs              = onemats->MOs;
    int nroao                = sysinfo->nroao;
    long long int nrofint    = sysinfo->nrofint;
    long long int* sortcount = twomats->sortcount;
    double* intval           = twomats->intval;
    unsigned short* intnums  = twomats->intnums;

    long long int prec_mem = (long long int) nroao*(long long int) nroao* (long long int) nroao* (long long int) nroao;
    std::cout << "Need " << prec_mem*sizeof(double) << " bytes (" << prec_mem*sizeof(double)/1024/1024 << " MB) for precalculation\n";

    postHF->prec_ints =  new double[prec_mem];

    //4-index trans
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;
    int i,j,k,l;
    long long int prec_count = 0;
    for(i = 0; i < nroao; i++) {
        for(j = 0; j < nroao; j++) {
            for(k = 0; k < nroao; k++) {
                for(l = 0; l <  nroao; l++) {
                    postHF->prec_ints[prec_count] = mo2int_op(i, j, k, l,nroao, MOs, nrofint, sortcount, intval, intnums,&std::cout);
                    prec_count++;
                }
            }
        }
        std::cout << i << "\t"<<std::flush;
        if((i+1)%10==0) std::cout << "\n"<<std::flush;
    }
}


/*******************************************************************************
*                                                                             *
* calc_mo2int_bf                                                              *
*                                                                             *
*                                                                             *
* bruit force claculation of 2el integral in MO-basis                         *
* (4 index transformation)  (chemist notation)                                *
*******************************************************************************/

//Inline functions for two electron integtral permutations

double perm_all(unsigned short a, unsigned short b,
                unsigned short c, unsigned short d, double integral,
                double* MOs, int off_i, int off_j, int off_k, int off_l){
    double inte = 0.;
    //J-Terme
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+c]*MOs[off_l+d]*integral; //1
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+d]*MOs[off_l+c]*integral; //2
    inte += MOs[off_i+b]*MOs[off_j+a]*MOs[off_k+c]*MOs[off_l+d]*integral; //3
    inte += MOs[off_i+b]*MOs[off_j+a]*MOs[off_k+d]*MOs[off_l+c]*integral; //4
    inte += MOs[off_i+c]*MOs[off_j+d]*MOs[off_k+a]*MOs[off_l+b]*integral; //5
    inte += MOs[off_i+d]*MOs[off_j+c]*MOs[off_k+a]*MOs[off_l+b]*integral; //6
    inte += MOs[off_i+c]*MOs[off_j+d]*MOs[off_k+b]*MOs[off_l+a]*integral; //7
    inte += MOs[off_i+d]*MOs[off_j+c]*MOs[off_k+b]*MOs[off_l+a]*integral; //8

    return(inte);
}

double perm_1234(unsigned short a, unsigned short b,
                 unsigned short c, unsigned short d, double integral,
                 double* MOs, int off_i, int off_j, int off_k, int off_l){
    double inte = 0.;
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+c]*MOs[off_l+d]*integral; //1
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+d]*MOs[off_l+c]*integral; //2
    inte += MOs[off_i+b]*MOs[off_j+a]*MOs[off_k+c]*MOs[off_l+d]*integral; //3
    inte += MOs[off_i+b]*MOs[off_j+a]*MOs[off_k+d]*MOs[off_l+c]*integral; //4

    return(inte);
}

double perm_1256(unsigned short a, unsigned short b,
                 unsigned short c, unsigned short d, double integral,
                 double* MOs, int off_i, int off_j, int off_k, int off_l){
    double inte = 0.;
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+c]*MOs[off_l+d]*integral; //1
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+d]*MOs[off_l+c]*integral; //2
    inte += MOs[off_i+c]*MOs[off_j+d]*MOs[off_k+a]*MOs[off_l+b]*integral; //5
    inte += MOs[off_i+d]*MOs[off_j+c]*MOs[off_k+a]*MOs[off_l+b]*integral; //6

    return(inte);
}

double perm_15(unsigned short a, unsigned short b,
               unsigned short c, unsigned short d, double integral,
               double* MOs, int off_i, int off_j, int off_k, int off_l){
    double inte = 0.;
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+c]*MOs[off_l+d]*integral; //1
    inte += MOs[off_i+c]*MOs[off_j+d]*MOs[off_k+a]*MOs[off_l+b]*integral; //5

    return(inte);
}

double perm_1(unsigned short a, unsigned short b,
              unsigned short c, unsigned short d, double integral,
              double* MOs, int off_i, int off_j, int off_k, int off_l){
    double inte = 0.;
    inte += MOs[off_i+a]*MOs[off_j+b]*MOs[off_k+c]*MOs[off_l+d]*integral; //1

    return(inte);
}


double calc_mo2int_bf(int i, int j, int k, int l, int nroao, double* MOs,
                      long long int nrofint, long long int* sortcount,  double* intval,
                      unsigned short* intnums){
    double inte = 0.;

    int off_i = nroao*i;
    int off_j = nroao*j;
    int off_k = nroao*k;
    int off_l = nroao*l;


    //PERM_1
    for(long long int x = 0; x < sortcount[0]; x++)
        inte += perm_1(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3],
                       intval[x], MOs, off_i, off_j, off_k, off_l);


    //PERM_15
    for(long long int x = sortcount[0]; x < sortcount[1]; x++)
        inte += perm_15(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3],
                        intval[x], MOs, off_i, off_j, off_k, off_l);

    //PERM_1234
    for(long long int x = sortcount[1]; x < sortcount[2]; x++)
        inte += perm_1234(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3],
                          intval[x], MOs, off_i, off_j, off_k, off_l);


    //PERM_1256
    for(long long int x = sortcount[2]; x < sortcount[3]; x++)
        inte += perm_1256(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3],
                          intval[x], MOs, off_i, off_j, off_k, off_l);


    //PERM_ALL
    for(long long int x = sortcount[3]; x < nrofint; x++)
        inte += perm_all(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3],
                         intval[x], MOs, off_i, off_j, off_k, off_l);


    return(inte);
}

/*******************************************************************************
*                                                                             *
* calc_mo2int_op                                                              *
*                                                                             *
*                                                                             *
* optimized     4 index transformation                                        *
* should be called with i being the slowest varying index                     *
* (4 index transformation)  (chemist notation)                                *
*******************************************************************************/


//TRANSFORM THE FIRST INDES AND BUILD ALL PERMUTATIONS OF THE REMAINING INDICES => AVOIDS K^4 MEMORY STORAGE
void trans_I(int i, double* i____Trans, int nroao, double* MOs,
             long long int nrofint, long long int* sortcount,  double* intval,
             unsigned short* intnums){

    //STATIC VARIABLES
    static long long int sM3 = (long long int) nroao *(long long int) nroao *(long long int) nroao;
    static long long int fac_d = nroao * nroao; //!!!!MAKE b the fast running index => will be transfromed next
    static long long int fac_c = nroao;

    //erase old I;
    for(long long int x = 0; x < sM3; x++) i____Trans[x] = 0.;

    int off_i = nroao*i;


    ///!!!!!!!!!!!DO TRANSFORMATION and PERMUTATIONS!!!!!!!!!!!!!!!

    //PERM_1
    for(long long int x = 0; x < sortcount[0]; x++) {
        //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
        i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
    }

    //PERM_15
    for(long long int x = sortcount[0]; x < sortcount[1]; x++) {
        //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
        i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
        i____Trans[intnums[x*4+3]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
    }

    //PERM_1234
    for(long long int x = sortcount[1]; x < sortcount[2]; x++) {
        //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
        i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
        i____Trans[intnums[x*4+1]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc)
        i____Trans[intnums[x*4+0]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(3)    (ba|cd)
        i____Trans[intnums[x*4+0]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(4)    (ba|dc)
    }


    //PERM_1256
    for(long long int x = sortcount[2]; x < sortcount[3]; x++) {
        //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
        i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
        i____Trans[intnums[x*4+1]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc)
        i____Trans[intnums[x*4+3]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
        i____Trans[intnums[x*4+2]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(6)    (dc|ab)
    }


    //PERM_ALL
    for(long long int x = sortcount[3]; x < nrofint; x++) {
        //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
        i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
        i____Trans[intnums[x*4+1]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc)
        i____Trans[intnums[x*4+0]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(3)    (ba|cd)
        i____Trans[intnums[x*4+0]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(4)    (ba|dc)
        i____Trans[intnums[x*4+3]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
        i____Trans[intnums[x*4+2]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(6)    (dc|ab)
        i____Trans[intnums[x*4+3]+intnums[x*4+1]*fac_c+intnums[x*4+0]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(7)    (cd|ba)
        i____Trans[intnums[x*4+2]+intnums[x*4+1]*fac_c+intnums[x*4+0]*fac_d] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(8)    (dc|ba)
    }
}



double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
                 long long int nrofint, long long int* sortcount,  double* intval,
                 unsigned short* intnums, std::ostream* outf){

    static double* dmem;      //pointer for memory allocation
    static double* i____Trans; //transformed ints for only i       (nroao^3)
    static double* ij___Trans; //transformed ints for only ij      (nroao^2)
    static double* ijk__Trans; //transfodmed ints for only ijk     (nroao  )

    static int alt_i = -1;
    static int alt_j = -1;
    static int alt_k = -1;

    static long long int stati = 0, statj=0, statk = 0;


    static long long int fac_d = nroao * nroao; //!!!!MAKE j the fast running index => will be transfromed next
    static long long int fac_c = nroao;


    if(dmem == NULL) {
        *outf << "...............................................................................\n";
        *outf << "Memory allocation for improved 4 index transformation\n";
        long long int M = nroao;
        long long int dM = M*M*M+M*M+M;
        *outf << "Need " << dM*8 << " bytes of Memory for one index  at a time transformation\n";
        outf->flush();
        dmem = new double[dM];
        long long int inc = 0;
        i____Trans = &(dmem[inc]); inc += M*M*M;
        ij___Trans = &(dmem[inc]); inc += M*M;
        ijk__Trans = &(dmem[inc]); inc += M;
        *outf << "Done \n";
        outf->flush();
        *outf << "...............................................................................\n";
    }

    //Print statistcs
    if(i == -1) {
        *outf << "\n...............................................................................\n";
        *outf << "mo2int_op transformation statistics: \n";
        *outf << "\nmo2int_op: # i trans so far " << stati << "\n";
        *outf <<   "mo2int_op: # j trans so far " << statj << "\n";
        *outf <<   "mo2int_op: # k trans so far " << statk << "\n";
        *outf << "...............................................................................\n";
        outf->flush();
        return(0);
    }



    if(i != alt_i) {
        //Make transformation for i
        trans_I( i,  i____Trans,  nroao,  MOs, nrofint, sortcount,   intval, intnums);
        alt_i = i;
        //Statistics
        stati++;
    }

    if(j != alt_j) {
        //Make transformation    j
        for(int d = 0; d < nroao; d++) {
            for(int c = 0; c < nroao; c++) {
                ij___Trans[d*nroao+c] = 0.;
                long long int ioff_b = fac_c*c +  fac_d*d;
                int off_j = nroao*j;
                for(int b = 0; b < nroao; b++)
                    ij___Trans[d*nroao+c] += i____Trans[ioff_b+b]*MOs[off_j+b];
            }
        }
        alt_j = j;
        //Statistics
        statj++;
    }

    if(k != alt_k) {
        //Make transformation for k
        for(int d = 0; d < nroao; d++) {
            int off_k = k*nroao;
            ijk__Trans[d] = 0.;
            for(int c = 0; c < nroao; c++)
                ijk__Trans[d] += ij___Trans[d*nroao+c]*MOs[off_k+c];
        }
        alt_k = k;
        //Statistics
        statk++;
    }


    //Make transfromation for l
    double inte = 0.;
    int off_l = nroao*l;
    for(int d = 0; d < nroao; d++)
        inte += ijk__Trans[d]*MOs[off_l+d];

    return(inte);
}
