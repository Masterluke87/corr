#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "motrans.h"


void build_FMo(int nroao,
               double* Fmat,
               double* FMo,
               double* MOs){

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



void read_transform_ri(std::string prefix,  //prefix to find the file
                       int nroe,            //nr of electrons
                       int nroao,           //nr of bsf in aobasis
                       int naux_2,          //nr of aux basis functions
                       double* MOs,         //orbitals for transformation
                       double* Bia)        //output - transformed array containing mointegrals
{

    int nocc = nroe/2;
    int nvir = nroao-nroe/2;
    double* BPQ = new double[naux_2*nroao*nroao];
    double* BiQ = new double[naux_2*nocc*nroao];

    memset(BiQ,0,sizeof(double)*naux_2*nocc*nroao);
    memset(Bia,0,sizeof(double)*naux_2*nocc*nvir);

    std::ifstream datf;
    datf.open(prefix+".rimp2");
    datf.read((char*) BPQ, naux_2*nroao*nroao*sizeof(double));
    datf.close();

    for (int Q=0; Q<naux_2; Q++)
        for (int i=0; i<nocc; i++)
            for (int mu=0; mu<nroao; mu++)
                for (int nu=0; nu<nroao; nu++)
                    BiQ[Q*nocc*nroao + i*nroao + nu] +=  MOs[i*nroao+mu]*BPQ[Q*nroao*nroao + mu*nroao + nu];
    for (int Q=0; Q<naux_2; Q++)
        for (int a=0; a<nvir; a++)
            for (int i =0; i<nocc; i++)
                for (int nu=0; nu<nroao; nu++)
                    Bia[Q*nocc*nvir+i*nvir+a]+= MOs[(a+nroe/2)*nroao + nu]*BiQ[Q*nocc*nroao + i*nroao + nu ];


    delete[] BiQ;
    delete[] BPQ;

}

void transform_ri(int nroe,            //nr of electrons
                  int nroao,           //nr of bsf in aobasis
                  int naux_2,          //nr of aux basis functions
                  double* MOs,         //orbitals for transformation
                  double* BPQ,          //input - calculated b^Q_pq
                  double* Bia)         //output - transformed array containing mointegrals)
{
    int nocc = nroe/2;
    int nvir = nroao-nroe/2;
    double* BiQ = new double[naux_2*nocc*nroao];
    memset(BiQ,0,sizeof(double)*naux_2*nocc*nroao);
    memset(Bia,0,sizeof(double)*naux_2*nocc*nvir);
    for (int Q=0; Q<naux_2; Q++)
        for (int i=0; i<nocc; i++)
            for (int mu=0; mu<nroao; mu++)
                for (int nu=0; nu<nroao; nu++)
                    BiQ[Q*nocc*nroao + i*nroao + nu] +=  MOs[i*nroao+mu]*BPQ[Q*nroao*nroao + mu*nroao + nu];
    for (int Q=0; Q<naux_2; Q++)
        for (int a=0; a<nvir; a++)
            for (int i =0; i<nocc; i++)
                for (int nu=0; nu<nroao; nu++)
                    Bia[Q*nocc*nvir+i*nvir+a]+= MOs[(a+nroe/2)*nroao + nu]*BiQ[Q*nocc*nroao + i*nroao + nu ];
    delete[] BiQ;
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
