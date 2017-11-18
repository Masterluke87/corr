#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "ops_rhf.h"

using namespace std;


//Functions
double calc_r_ab(int a, int b, double* coord);
double calc_ion_rep(int nroa, double* coord, double* charges);
void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
void calc_mu_core(int nroa, double* coord, double* charges, double* point,
		  double* mu_core);
double calc_S12(int nroao, double* Smat, double* Som12, double* tmpmat,
		double* tmpvecs, double* tmpvals);
double calc_Sp12(int nroao, double* Smat, double* Sop12, double* tmpmat,
		 double* tmpvecs, double* tmpvals);
void   diag_Fmat(int nroao, double* Fmat, double* MOs,
		 double* MOens, double* Som12, double* tmpmat);
void   build_Fmat(int nroao, double* Fmat, double* Pmat, double* Hmat,
		  double* intval, unsigned short* intnums,
		  long long int* sortcount, long long int nrofint);
double Calc_e_el(int nroao, double* Fmat, double* Pmat, double* Hmat);
double calc_op_1el(int nroao, double* opmat, double* Pmat);

//Extern Functions
extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);
extern void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
extern void transform_MOs(int nroao, double *MOs, double* tmat, double* tmpvec);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_r_ab                                                       */
/*                                                                               */
/* caluculates the ion core coulomb repulsion                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_r_ab(int a, int b, double* coord){
  double r = sqrt(pow(coord[3*a+0]-coord[3*b+0],2)+
		  pow(coord[3*a+1]-coord[3*b+1],2)+
		  pow(coord[3*a+2]-coord[3*b+2],2));
  return(r);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_ion_rep                                                    */
/*                                                                               */
/* calculates the ion core coulomb repulsion                                     */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_ion_rep(int nroa, double* coord, double* charges){
  double ion_rep = 0.;
  for(int x = 0; x < nroa; x++){
    for(int y = x+1; y < nroa; y++){
      ion_rep += charges[x]*charges[y]/calc_r_ab(x,y,coord);
    }
  }

  return(ion_rep);
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_center_of_mass                                             */
/*                                                                               */
/* calculates the ion core center of mass                                         */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass){
  for(int x = 0; x < 3; x++) center_of_mass[x] = 0.;

  double mass_tot = 0.;

  for(int x = 0; x < nroa; x++){
    mass_tot += mass[x];
    for(int i = 0; i < 3; i++)
      center_of_mass[i] += mass[x] * coord[3*x+i];
  }

  for(int i = 0; i < 3; i++)
    center_of_mass[i] /= mass_tot;

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_mu_core                                                    */
/*                                                                               */
/* calculates the ion core dipole moment  at point                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_mu_core(int nroa, double* coord, double* charges, double* point,
		  double* mu_core){
  for(int x = 0; x < 3; x++) mu_core[x] = 0.;


  for(int x = 0; x < nroa; x++){
    for(int i = 0; i < 3; i++)
      mu_core[i] += charges[x] * (coord[3*x+i] - point[i]);
  }

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_S12                                                        */
/*                                                                               */
/* calculate S^-1/2 for symmetric orthogonalisation                              */
/* returns mimum eigenvalue of S                                                 */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_S12(int nroao, double* Smat, double* Som12, double* tmpmat,
		double* tmpvecs, double* tmpvals){

  diag_mat(nroao, Smat, tmpvals, tmpvecs);

  double min_val =10.;
  for(int x = 0; x < nroao; x++){
    if(min_val > tmpvals[x]) min_val = tmpvals[x];
  }

  for(int x = 0; x <  nroao*nroao; x++){
    Som12[x] = 0.; tmpmat[x] = 0.;
  }
  for(int x = 0; x <  nroao; x++)
    Som12[x*nroao+x] = 1./sqrt(tmpvals[x]);
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      for(int z = 0; z < nroao; z++)  //Adjungiert !!!!!!
        tmpmat[x*nroao+y] +=  Som12[z*nroao+y] * tmpvecs[z*nroao+x];
    }
  }
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      Som12[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
        Som12[x*nroao+y] += tmpvecs[z*nroao+y] * tmpmat[x*nroao+z];
    }
  }

  return(min_val);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_Sp12                                                       */
/*                                                                               */
/* calculate S^+1/2 for symmetric orthogonalisation of P                         */
/* returns mimum eigenvalue of S                                                 */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_Sp12(int nroao, double* Smat, double* Sop12, double* tmpmat,
		double* tmpvecs, double* tmpvals){

  diag_mat(nroao, Smat, tmpvals, tmpvecs);

  double min_val =10.;
  for(int x = 0; x < nroao; x++){
    if(min_val > tmpvals[x]) min_val = tmpvals[x];
  }

  for(int x = 0; x <  nroao*nroao; x++){
    Sop12[x] = 0.; tmpmat[x] = 0.;
  }
  for(int x = 0; x <  nroao; x++)
    Sop12[x*nroao+x] = sqrt(tmpvals[x]);
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      for(int z = 0; z < nroao; z++)  //Adjungiert !!!!!!
        tmpmat[x*nroao+y] +=  Sop12[z*nroao+y] * tmpvecs[z*nroao+x];
    }
  }
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      Sop12[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
        Sop12[x*nroao+y] += tmpvecs[z*nroao+y] * tmpmat[x*nroao+z];
    }
  }

  return(min_val);
}

/*******************************************************************************
 * This routine diagonalises the fock operator                                 *
 *                                                                             *
 ******************************************************************************/

void diag_Fmat(int nroao, double* Fmat, double* MOs,
	       double* MOens, double* Som12, double* tmpmat){
  symmortho_mat(nroao,Fmat,Som12,tmpmat);
  diag_mat(nroao,Fmat,MOens,MOs);
  transform_MOs(nroao, MOs, Som12,tmpmat);
}


/*******************************************************************************
 * This routine builds the new Pmat in a damped SCF                            *
 *                                                                             *
 ******************************************************************************/

double build_Pmat_dscf(int nroao, int nroe, double* Pmat, double* Pmat_old,
		       double* MOs, double damp){
  double max_diff = 0.;

  for(int x = 0; x < nroao*nroao; x++) Pmat[x] = 0.;

  for(int e = 0; e < nroe/2; e++){
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao; y++){
        Pmat[x*nroao+y] += 2.*MOs[e*nroao+y]*MOs[e*nroao+x];
      }
    }
  }

  for(int x = 0; x < nroao*nroao; x++){
    if(fabs(Pmat_old[x] - Pmat[x]) > max_diff) max_diff = fabs(Pmat_old[x] - Pmat[x]);
    Pmat[x] = (1.-damp)*Pmat[x] + damp*Pmat_old[x];
    Pmat_old[x] = Pmat[x];
  }

  return(max_diff);
}

/*******************************************************************************
 * Permutations for fock matrix build                                          *
 *                                                                             *
 ******************************************************************************/


inline void perm_all(unsigned short a, unsigned short b,
		     unsigned short c, unsigned short d, double integral,
                     double* Fmat, double* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
  Fmat[b*nroao+a] += Pmat[c*nroao+d]*integral; //3
  Fmat[b*nroao+a] += Pmat[d*nroao+c]*integral; //4
  Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5
  Fmat[d*nroao+c] += Pmat[a*nroao+b]*integral; //6
  Fmat[c*nroao+d] += Pmat[b*nroao+a]*integral; //7
  Fmat[d*nroao+c] += Pmat[b*nroao+a]*integral; //8

  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
  Fmat[b*nroao+c] -= 0.5*Pmat[a*nroao+d]*integral; //3
  Fmat[b*nroao+d] -= 0.5*Pmat[a*nroao+c]*integral; //4
  Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
  Fmat[d*nroao+a] -= 0.5*Pmat[c*nroao+b]*integral; //6
  Fmat[c*nroao+b] -= 0.5*Pmat[d*nroao+a]*integral; //7
  Fmat[d*nroao+b] -= 0.5*Pmat[c*nroao+a]*integral; //8
}

inline void perm_1234(unsigned short a, unsigned short b,
		      unsigned short c, unsigned short d, double integral,
                      double* Fmat, double* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
  Fmat[b*nroao+a] += Pmat[c*nroao+d]*integral; //3
  Fmat[b*nroao+a] += Pmat[d*nroao+c]*integral; //4

  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
  Fmat[b*nroao+c] -= 0.5*Pmat[a*nroao+d]*integral; //3
  Fmat[b*nroao+d] -= 0.5*Pmat[a*nroao+c]*integral; //4
}


inline void perm_1256(unsigned short a, unsigned short b,
		      unsigned short c, unsigned short d, double integral,
                      double* Fmat, double* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
  Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5
  Fmat[d*nroao+c] += Pmat[a*nroao+b]*integral; //6

  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
  Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
  Fmat[d*nroao+a] -= 0.5*Pmat[c*nroao+b]*integral; //6
}

inline void perm_15(unsigned short a, unsigned short b,
		    unsigned short c, unsigned short d, double integral,
                    double* Fmat, double* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5

  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
}

inline void perm_1(unsigned short a, unsigned short b,
		   unsigned short c, unsigned short d, double integral,
                   double* Fmat, double* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1

  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
}

/*******************************************************************************
 * Fock matrix build                                                           *
 *                                                                             *
 ******************************************************************************/
void build_Fmat(int nroao, double* Fmat, double* Pmat, double* Hmat,double* intval, unsigned short* intnums, long long int* sortcount, long long int nrofint){
  for(int mu = 0; mu < nroao; mu++){
    for(int nu = 0; nu < nroao;nu++){
      Fmat[mu*nroao+nu] = Hmat[mu*nroao+nu];
    }
  }


  //PERM_1
  for(long long int x = 0; x < sortcount[0]; x++)
    perm_1(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);


  //PERM_15
  for(long long int x = sortcount[0]; x < sortcount[1]; x++)
    perm_15(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);


  //PERM_1234
  for(long long int x = sortcount[1]; x < sortcount[2]; x++)
    perm_1234(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);


  //PERM_1256
  for(long long int x = sortcount[2]; x < sortcount[3]; x++)
    perm_1256(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);


  //PERM_ALL
  for(long long int x = sortcount[3]; x < nrofint; x++)
    perm_all(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);
}


/*******************************************************************************
 *  Calculation of electronic HF-energy                                        *
 *                                                                             *
 ******************************************************************************/

double Calc_e_el(int nroao, double* Fmat, double* Pmat, double* Hmat){
  double  E_tot_el = 0.;
  for(int x = 0; x < nroao;x++){
    for(int y = 0; y < nroao;y++){
      E_tot_el += Pmat[x*nroao+y]*(Fmat[y*nroao+x] + Hmat[y*nroao+x]);
    }
  }
  return(E_tot_el/2.);
}
/*******************************************************************************
 * calculation of one electrion expc. values                                   *
 *                                                                             *
 ******************************************************************************/

double   calc_op_1el(int nroao, double* opmat, double* Pmat){
  double op = 0.;
  for(int x = 0; x < nroao;x++){
    for(int y = 0; y < nroao;y++){
      op += Pmat[x*nroao+y]*opmat[y*nroao+x];
    }
  }
  return(op);
}
