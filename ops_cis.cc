/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_cis.cc                                                             *
 *                                                                              *
 * contains functions for especially needed for cis calculations                *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 ********************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "ops_cis.h"



//Extern Functions
extern double   calc_op_1el(int nroao, double* opmat, double* Pmat);

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
  static long long int fac_d = nroao * nroao;   //!!!!MAKE b the fast running index => will be transfromed next
  static long long int fac_c = nroao;

  //erase old I;
  for(long long int x = 0; x < sM3; x++) i____Trans[x] = 0.;

  int off_i = nroao*i;


  ///!!!!!!!!!!!DO TRANSFORMATION and PERMUTATIONS!!!!!!!!!!!!!!!

  //PERM_1
  for(long long int x = 0; x < sortcount[0]; x++){
    //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
    i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
  }

  //PERM_15
  for(long long int x = sortcount[0]; x < sortcount[1]; x++){
    //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!)
    i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
    i____Trans[intnums[x*4+3]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
  }

  //PERM_1234
  for(long long int x = sortcount[1]; x < sortcount[2]; x++){
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
  for(long long int x = sortcount[3]; x < nrofint; x++){
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

  static double* dmem;        //pointer for memory allocation
  static double* i____Trans;  //transformed ints for only i       (nroao^3)
  static double* ij___Trans;  //transformed ints for only ij      (nroao^2)
  static double* ijk__Trans;  //transfodmed ints for only ijk     (nroao  )

  static int alt_i = -1;
  static int alt_j = -1;
  static int alt_k = -1;

  static long long int stati = 0, statj=0, statk = 0;


  static long long int fac_d = nroao * nroao;   //!!!!MAKE j the fast running index => will be transfromed next
  static long long int fac_c = nroao;


  if(dmem == NULL){
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
  if(i == -1){
    *outf << "\n...............................................................................\n";
    *outf << "mo2int_op transformation statistics: \n";
    *outf << "\nmo2int_op: # i trans so far " << stati << "\n";
    *outf <<   "mo2int_op: # j trans so far " << statj << "\n";
    *outf <<   "mo2int_op: # k trans so far " << statk << "\n";
    *outf << "...............................................................................\n";
    outf->flush();
    return(0);
  }



  if(i != alt_i){
    //Make transformation for i
    trans_I( i,  i____Trans,  nroao,  MOs, nrofint, sortcount,   intval, intnums);
    alt_i = i;
    //Statistics
    stati++;
  }

  if(j != alt_j){
    //Make transformation    j
    for(int d = 0; d < nroao; d++){
      for(int c = 0; c < nroao; c++){
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

  if(k != alt_k){
    //Make transformation for k
    for(int d = 0; d < nroao; d++){
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


/*******************************************************************************
 *                                                                             *
 * calc_1p_op_cis_od                                                           *
 *                                                                             *
 *                                                                             *
 * one particle operator between cis basis functions   (off doagonal elements) *
 *******************************************************************************/

double calc_1p_op_cis_od(int a, int b, double* MOs, int nroao, double* mat,
			 double* tmpvec){
  double op = 0.;

  for(int x = 0; x < nroao; x++){
    tmpvec[x] = 0.;
    for(int y = 0; y < nroao; y++)
      tmpvec[x] += mat[x*nroao+y] * MOs[b*nroao+y];
    op += tmpvec[x]*MOs[a*nroao+x];
  }

  return(op);
}


/*******************************************************************************
 *                                                                             *
 * calc_1p_op_cis_d                                                            *
 *                                                                             *
 *                                                                             *
 * one particle operator between cis basis functions  (diagonal elements)      *
 *******************************************************************************/

double calc_1p_op_cis_d(int i, int f, double* MOs, int nroao, double* opmat,
			 double* Pmat, int nroe){
  for(int x = 0; x <nroao*nroao; x++) Pmat[x] = 0.;

  for(int e = 0; e < nroe/2; e++){
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao; y++)
	Pmat[x*nroao+y] += 2.*MOs[e*nroao+x]*MOs[e*nroao+y];
    }
  }

  for(int x = 0; x < nroao; x++){
    for(int y = 0; y < nroao; y++)
      Pmat[x*nroao+y] +=  -MOs[i*nroao+x]*MOs[i*nroao+y]+MOs[f*nroao+x]*MOs[f*nroao+y];
  }

  double opval = calc_op_1el( nroao,  opmat,  Pmat);
  return(opval);
}


/*******************************************************************************
 *                                                                             *
 * calc_mu_mat_cis                                                             *
 *                                                                             *
 *                                                                             *
 * calculates dipole matrix in cis_basis                                       *
 *******************************************************************************/


void calc_mu_mat_cis(int cis_size, int nroao, int omo, int umo, int llim, int nroe,
		     double *cismat, double* mumat, double mucore,
		     double* MOs, double* Pmat, double* tmpvec_ao){


  long long int Cis_size = cis_size;
  for(long long int x = 0; x < Cis_size*Cis_size; x++)
    cismat[x] = 0.;


  cismat[0] = -calc_1p_op_cis_d( 0, 0,  MOs,  nroao,  mumat, Pmat, nroe)+mucore;

  for(long long int x = 1; x < Cis_size; x++){
    long long int i = (x-1)/umo+llim;
    long long int f = (x-1)%umo+omo+llim;
    cismat[x] =  sqrt(2.)*calc_1p_op_cis_od( i, f, MOs,  nroao, mumat, tmpvec_ao);
  }

  for(long long int x = 1 ; x < Cis_size; x++){
    long long int i1 = (x-1)/umo+llim;
    long long int f1 = (x-1)%umo+omo+llim;
    for(long long int y = x ; y < Cis_size; y++){
      long long int i2 = (y-1)/umo+llim;
      long long int f2 = (y-1)%umo+omo+llim;
      if(i1==i2 && f1!=f2)
        cismat[x*Cis_size+y] =  -calc_1p_op_cis_od( f2, f1, MOs,  nroao, mumat, tmpvec_ao);
      if(i1!=i2 && f1==f2)
        cismat[x*Cis_size+y] =   calc_1p_op_cis_od( i1, i2, MOs,  nroao, mumat, tmpvec_ao);
      if(y == x){
	cismat[x*Cis_size+y] =  -calc_1p_op_cis_d( i1, f1,  MOs,  nroao,  mumat, Pmat, nroe)+mucore;
      }
    }
  }

  for(long long int x = 0; x < Cis_size; x++){
    for(long long int y = x; y < Cis_size; y++){
      cismat[y*Cis_size+x] = cismat[x*Cis_size+y];
    }
  }
}

/*******************************************************************************
 * ROUTINES TO MAP i j k l for needed integrals to 2D fiel indecies x,y        *
 *                                                                             *
 * Needed integrals (CIS, MP2, CIS(D))  (o = occ. MO, v = virt. MO)            *
 *                                                                             *
 * (oo|vv), (ov|ov), (ov|vv), (oo|ov)                                          *
 *                                                                             *
 *   => We need (i j | k l )   i = llim...homo                                 *
 *                             j = i   ...homo                                 *
 *                             k = llim...ulim                                 *
 *                             l = lumo or k ...ulim                           *
 *                                                                             *
 *  We will use to indices x(i,j) and y(k,l)                                   *
 *  with conditions i<=j, k<=l, precalculated integrals are stored             *
 *  in intmo_pc[x*y_max+y] with y_max = y(ulim,ulim)                           *
 *                                                                             *
 *     x(i,j) is if i-llim == 0 ->  x(i,j) = j-llim                            *
 *                                                                             *
 *               else           ->  x(i,j) = Prod_m=0^m<i-llim (omo+umo-m)     *
 *                                         +j-llim                             *
 *                                                                             *
 *    y(k,l) is  k <= homo      ->  y(k,l) = (k-llim)*umo+l-llim-omo           *
 *                                                                             *
 *               else           ->  y(k,l) = omo*umo+                          *
 *                                           SUM_m=0^m<k-llim-omo (umo-m)      *
 *                                           +l-k                              *
 *                                                                             *
 *                                                                             *
 *   ALLWAYS BRING i j k l into the right order (i<=j, k<=l) before calling    *
 *   one of the following routines!!!!!!!!                                     *
 *                                                                             *
 *                                                 Tillmann Klamroth 2004      *
 *******************************************************************************/

//x(i,j)
int calc_mo2el_ind_x(int i, int j, int omo, int umo, int llim){
  // if( i-llim == 0) return(j-llim);
  int x = 0;
  for(int m = 0; m <  i-llim;m++) x += omo+umo-m;
  return(x+j-i);
}


//y(k,l)
int calc_mo2el_ind_y(int k, int l, int omo, int umo, int llim){
  if(k-llim <= omo) return((k-llim)*umo+l-llim-omo);
  int y = 0;
  for(int m = 0; m < k-llim-omo;m++) y += (umo-m);
  return(y+omo*umo+l-k);
  //  return(y+(k-llim)*umo+l-llim-omo);
}


/*******************************************************************************
 *                                                                             *
 * precalc_ints_sd                                                             *
 *                                                                             *
 *                                                                             *
 * repcalculates all integrals needed for CIS, MP2, CIS(D)                     *
 * returns pointer (double) to sorted integrals                                *
 *******************************************************************************/


double* precalc_ints_sd(int omo, int umo, int llim,
			int nroao, double* MOs,
			long long int nrofint, long long int* sortcount,  double* intval,
			unsigned short* intnums, std::ofstream* outf){
  *outf << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  *outf << "precalc_ints_sd:\n";

  //Get needed Memory
  long long int x_max = calc_mo2el_ind_x(omo-1, omo+umo, omo, umo, 0);
  long long int y_max=  calc_mo2el_ind_y(omo+umo,omo+umo, omo,  umo, 0);
  *outf << "Needed storage is " << x_max << " x " << y_max << " => " << x_max*y_max*8 << " bytes \n";

  double* prec_ints = new double[x_max*y_max];
  outf->flush();
  *outf<< "Loop over O: for (Oo|ov), (Oo|vv), (Ov|ov), (Ov|vv): \n";

  long long int prec_count = 0;

  for(int i = llim; i < llim+omo; i++){
    for(int j = i; j < llim+omo+umo; j++){
      for(int k = llim; k < llim+omo+umo; k++){
        int K = k;
        if(K < llim+omo) K = llim+omo;
        for(int l = K; l < llim+omo+umo;l++){
	  prec_ints[prec_count++] = mo2int_op(i, j, k, l,
					      nroao, MOs, nrofint, sortcount, intval, intnums, outf);
	}
      }
    }
    *outf << i << "\t";
    if((i+1)%10==0) *outf << "\n";
    outf->flush();
  }
  *outf << "\nInts claculated: " << prec_count << "\n";
  *outf << "\nDone\n";
  outf->flush();
  //Print out stat
  mo2int_op(-1, -1, -1, -1,
	    nroao, MOs, nrofint, sortcount, intval, intnums, outf);

  return(prec_ints);
}

/*******************************************************************************
 *                                                                             *
 * get_precalc_ints_sd                                                         *
 *                                                                             *
 * extract precalculated integrals  safe                                       *
 *******************************************************************************/


double get_precalc_ints_sd(int i, int j, int k, int l,
			   int omo, int umo, int llim, double* prec_ints){
  int tmp;
  //REORDER

  if(i>j){tmp = i; i = j; j = tmp;}
  if(k>l){tmp = k; k = l; l = tmp;}

  if(i>=llim+omo)
    { tmp = i; i = k; k = tmp;
      tmp = j; j = l; l = tmp; }


  int x_max  = calc_mo2el_ind_x(omo, omo+umo-1, omo,umo, 0);
  int y_max  = calc_mo2el_ind_y(omo+umo,omo+umo, omo,  umo, 0);

  int x = calc_mo2el_ind_x(i, j, omo, umo, llim);
  int y = calc_mo2el_ind_y(k, l, omo, umo, llim);
  //BOUNCE CHECK
  if(x > x_max || y > y_max) {
    std::cerr << "Missing integral (" << i << " " << j << "|" << k << " " << l << ")\n";
    exit(0);
  }

  return(prec_ints[x*y_max+y]);
}

/*******************************************************************************
 *                                                                             *
 * get_precalc_ints_sd                                                         *
 *                                                                             *
 * extract precalculated integrals fast ovov type                              *
 * special routine for RO_CIS(D)
 *******************************************************************************/

int* init_prex( int omo, int umo, int llim){
  int*  prex = new int[omo];
  for(int x = 0; x < omo; x++)
    prex[x] = calc_mo2el_ind_x(x+llim, 0, omo, umo, llim);
  return(prex);
}

double get_precalc_ints_ovov(int i, int j, int k, int l,
			     int omo, int umo, int llim, double* prec_ints){
  //ONLY ONE TIME INIT
  static  int y_max  = calc_mo2el_ind_y(omo+umo,omo+umo, omo,  umo, 0);
  static  int* prex  =  init_prex( omo,umo,llim);

  int x = prex[i-llim]+j;         //calc_mo2el_ind_x(i, j, omo, umo, llim);
  int y = (k-llim)*umo+l-llim-omo;

  return(prec_ints[x*y_max+y]);
}
