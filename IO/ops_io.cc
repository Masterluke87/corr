#include <fstream>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "ops_io.h"

using namespace std;

//Extern Functions




/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      STATUS                                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void status(ofstream* outf){
  char str1[1024];

  size_t len = 1024;
  time_t curr_time;
  time(&curr_time);

  gethostname(str1, len);

  *outf << "Host: " << str1 << ", Date: " << ctime(&curr_time);
  outf->flush();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      REM COM                                                  */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int rem_com(char* filename, char* streamstring, int string_length){
  const char com_B = '#';
  const char com_E = '\n';

  int pos = 0;
  char cc;

  ifstream inf(filename);

  while(inf.get(cc)&& pos < string_length-1){
    if(cc != com_B)
      streamstring[pos++] = cc;
    else{
      while(cc != com_E && inf.get(cc));
      streamstring[pos++] = com_E;
    }
  }
  streamstring[pos] = 0;
  if(pos == string_length-1){
    cerr << "Buffer size exceeded !\n"; exit(0);
  }
  return(strlen(streamstring));
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                get_sys_size(char* sysfile, int* nroao, int* nroa,             */
/*                              long long int* nrofint)                          */
/*                                                                               */
/* Read the system sizes form sysfile.                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void get_sys_size(std::string sysfile, int* nroao, int* nroa, long long int* nrofint){

  ifstream inf(sysfile);

  inf.read((char *) nroao,   sizeof(int));
  inf.read((char *) nroa,    sizeof(int));
  inf.read((char *) nrofint, sizeof(long long int));

  inf.close();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_sys(........)                                             */
/*                                                                               */
/*                                                                               */
/* Read the system data                                                          */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_sys(std::string sysfile, double* coord, double* charges, double* mass,
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy,
	      double *Dz, long long int* sortcount, double* intval,
	      unsigned short* intnums){
  ifstream datf(sysfile);

  int nroao, nroa;
  long long int nrofint;

  //SYSTEM DATA
  datf.read((char *) &nroao , sizeof(int));
  datf.read((char *) &nroa  , sizeof(int));
  datf.read((char *) &nrofint,  sizeof(long long int));
  datf.read((char *) coord  , sizeof(double)*3*nroa);
  datf.read((char *) charges, sizeof(double)*nroa);
  datf.read((char *) mass,    sizeof(double)*nroa);

  //ONEL EL INTEGRAL DATA
  datf.read((char * ) Hmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Tmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Smat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dx    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dy    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dz    , sizeof(double)*nroao*nroao);

  //TWO EL INTEGRAL DATA
  datf.read((char *) sortcount, sizeof(long long int)*4);
  datf.read((char *) intval,    sizeof(double)*nrofint);
  datf.read((char *) intnums,   sizeof(unsigned short)*nrofint*4);

  datf.close();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_sys_1el(........)                                         */
/*                                                                               */
/*                                                                               */
/* Read the system data, only 1 electron integrals                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_sys_1el(char* sysfile, double* coord, double* charges, double* mass,
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy,
	      double *Dz){
  ifstream datf(sysfile);

  int nroao, nroa;
  long long int nrofint;

  //SYSTEM DATA
  datf.read((char *) &nroao , sizeof(int));
  datf.read((char *) &nroa  , sizeof(int));
  datf.read((char *) &nrofint,  sizeof(long long int));
  datf.read((char *) coord  , sizeof(double)*3*nroa);
  datf.read((char *) charges, sizeof(double)*nroa);
  datf.read((char *) mass,    sizeof(double)*nroa);

  //ONEL EL INTEGRAL DATA
  datf.read((char * ) Hmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Tmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Smat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dx    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dy    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dz    , sizeof(double)*nroao*nroao);

  datf.close();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               write_wav_HF(........)                                          */
/*                                                                               */
/*                                                                               */
/* write  HF wave_function                                                       */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void write_wav_HF(std::string wavfile, int nroao, double* MOens, double* MOs){
  ofstream outf(wavfile);

  outf.write((char *) &nroao, sizeof(int));
  outf.write((char *) MOens,  sizeof(double)*nroao);
  outf.write((char *) MOs,    sizeof(double)*nroao*nroao);

  outf.close();

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_wav_HF(........)                                          */
/*                                                                               */
/*                                                                               */
/* Read  HF wave_function                                                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void read_wav_HF(std::string wavfile, int nroao, double* MOens, double* MOs){
  ifstream inf(wavfile);
  int real_nroao;

  inf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao != nroao){
    cerr << "Wrong HF wavefunction size in read_wav_HF!\n"; exit(3);
  }

  //DEBUG 04052012TK
  //cerr << real_nroao << " " << nroao << "\n";

  inf.read((char *) MOens,  sizeof(double)*nroao);
  inf.read((char *) MOs,    sizeof(double)*nroao*nroao);

  inf.close();
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                output_matrix(...)                                             */
/*                                                                               */
/*                                                                               */
/* Output lower triangular of symmetric matrix                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void   output_matrix(double* mat,int nr_of_col, int nrop, ofstream *outf){
  char dumchar[64];
  int nr_of_blocks = (int) ceil((double)nrop/(double)nr_of_col);
  for(int l = 0; l <  nr_of_blocks; l++){
    int upper_bound = (l+1)*nr_of_col;
    if(upper_bound > nrop) upper_bound = nrop;
    for(int x = l*nr_of_col; x < upper_bound; x++)
      *outf << "\t\t" << x+1 ;
    *outf << "\n\n";
    //loops over centers
    for(int y = l*nr_of_col; y < nrop; y++){
      *outf << y+1 << "\t" ;
      for(int x = l*nr_of_col; x < upper_bound && x <= y; x++){
      sprintf(dumchar,"%+.6e",mat[x*nrop+y]);
           *outf << dumchar << "\t";
      }
      *outf << "\n";
    }
    *outf << "\n";
  }
  outf->flush();
}
