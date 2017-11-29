#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <complex>
#include <string.h>
#include "ops_mat.h"




//Functions

//Extern Functions
extern "C" void  dsyev_(char* JOBZ, char*  UPLO,int* N, double* A, int* LDA,
                        double* W, double* WORK, int* LWORK, int*  INFO );

/*******************************************************************************
 * Matrix diagonalization                                                      *
 *                                                                             *
 ******************************************************************************/

void   diag_mat(int nroao, double* mat, double* vals, double* vecs){
  static int N = -1 ;
  static double* WORK;

  int LDA = nroao, LWORK=3*nroao-1, INFO;

  if(N != nroao){
    if(WORK != NULL){
      delete [] WORK;
     }

    WORK = new double[LWORK];
    N = nroao;
  }

  char JOBZ, UPLO;
  double* A    = vecs;
  double* W    = vals;
  JOBZ = 'V';
  UPLO = 'U';

  for(int x = 0; x <  nroao*nroao; x++)
    A[x] = mat[x];
  dsyev_( &JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO );
}

/*******************************************************************************
 * matrix transformation of mat with tmat (Som12)                              *
 *                                                                             *
 *  matneu = tmat mat tmat                                                     *
 *******************************************************************************/

void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat){
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      dummat[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
	dummat[x*nroao+y] += tmat[x*nroao+z] * mat[z*nroao+y];
    }
  }
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      mat[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
        mat[x*nroao+y] += tmat[z*nroao+y] * dummat[x*nroao+z];
    }
  }
}

/*******************************************************************************
 * Transformation of  MO's back to AO-Basis                                    *
 *                                                                             *
 ******************************************************************************/

void transform_MOs(int nroao, double *MOs, double* tmat, double* tmpvec){

  for(int x = 0; x < nroao;x++){
    for(int y = 0; y < nroao; y++){
      tmpvec[y] = 0.;
      for(int z = 0; z < nroao; z++){
        tmpvec[y] += tmat[y*nroao+z]*MOs[x*nroao+z];
      }
    }
    for(int y = 0; y < nroao; y++)
      MOs[x*nroao+y] = tmpvec[y];
  }
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     SIMPLE MATRIX ROUTINES                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f){
  int x,y,z;
  #pragma omp parallel for private(y,z)
  for(x = 0; x < np; x++){
    for(y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
    }
  }
}


void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f){
  int x,y,z;
  #pragma omp parallel for private(y,z)
  for(x = 0; x < np; x++){
    for(y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
    }
  }
}



void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
  long long int x,y,z;
  #pragma omp parallel for private(y,z)
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
    }
  }
}


void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
  long long int x,y,z;
  #pragma omp parallel for private(y,z)
  for(x = 0; x < np; x++){
    for(y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
    }
  }
}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     BASIS TRANSFORMAITIONS                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// FOR ALL ROUTINES
// DIR = 1 means:   mat = trans * mat * trans^T
// else             mat = trans^T*mat*trans

/*
void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir){

  if(dir == 0) //CALC mat*trans
    mat_mat( np,  mat,trans, tmpmat);
  else         //CALC trans*mat
    mat_mat( np, trans, mat, tmpmat);


  if(dir == 0) //CALC trans^T*(mat*trans)
    mat_T_mat( np, trans, tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    mat_mat_T( np, tmpmat, trans, mat);

}
*/
// For large matrices
/*
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir){

  if(dir == 0) //CALC mat*trans
    mat_mat( np,  mat,trans, tmpmat);
  else         //CALC trans*mat
    mat_mat( np, trans, mat, tmpmat);


  if(dir == 0) //CALC trans^T*(mat*trans)
    mat_T_mat( np, trans, tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    mat_mat_T( np, tmpmat, trans, mat);

}
*/

void pmv(double* mat, double* vi, double* vo, int nroao){
  for(int x = 0; x < nroao; x++){
    vo[x] = 0.;
    for(int y = 0; y < nroao; y++)
      vo[x] += mat[x*nroao+y]*vi[y];
  }
}

