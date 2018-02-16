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
#include <iomanip>
#include "ops_io.h"


//Extern Functions




/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      STATUS                                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void print_header(int* nroe,int* nroa,int* nroao,int* naux_1,int* naux_2,long long int* nrofint,
                  long long int* nrofaux,long long int* nrofaux2,
                  double* coord,double* charges,double* zeff,
                  std::string* basisNameOB,std::string* basisNameJK,std::string* basisNameRI){
  char str1[1024];

  size_t len = 1024;
  time_t curr_time;
  time(&curr_time);

  gethostname(str1, len);
  std::cout << "++++++++++++++++++++++++++++++++++++++\n";
  std::cout << "+A PSI4/MPI/OPENMP/CUDA - MP2 PROGRAM+\n";
  std::cout << "++++++++++++++++++++++++++++++++++++++\n";

  std::cout << "Host: " << str1 << "\nDate: " << ctime(&curr_time);
  std::cout << "\nSYSTEMDATA" << '\n';
  std::cout << "----------" << "\n\n";
  std::cout<<std::setw(-2)<<"Z" <<std::setw(10)<<"zeff"<<std::setw(10)<<"x" <<std::setw(10) << "y" <<std::setw(10)<< "z" <<'\n';
  for (size_t i = 0; i < *nroa; i++) {
        std::cout<<std::setw(-2)<<charges[i]<<std::setw(10)<<zeff[i]<<std::setw(10)<< coord[i*3+0]<<std::setw(10) << coord[i*3+1]<<std::setw(10)<< coord[i*3+2]<<'\n';
  }

  std::cout << "nroe          :" << *nroe<<'\n';
  std::cout << "nroa          :" << *nroa<<'\n';
  std::cout << "nroao         :" << *nroao<<'\n';
  std::cout << "naux-JK       :" << *naux_1<<'\n';
  std::cout << "naux-RI       :" << *naux_2<<'\n';
  std::cout << "nrofint       :" << *nrofint<<'\n';
  std::cout << "nrofaux(JK)   :" << *nrofaux<<'\n';
  std::cout << "nrofaux(RI)   :" << *nrofaux2<< '\n';
  std::cout << "basisName(OB) :" << *basisNameOB<<'\n';
  std::cout << "basisName(JK) :" << *basisNameJK<<'\n';
  std::cout << "basisName(RI) :" << *basisNameRI<<'\n';

  std::cout.flush();
}

void read_system(std::string filename,int* nroe,int* nroa,
                  int* nroao, int* naux_1, int* naux_2,
                  long long int* nrofint,long long int* nrofaux,long long int* nrofaux2,
                  double** coord,double** charges,double** zeff, double** mass,
                  std::string* basisNameOB, std::string* basisNameJK, std::string* basisNameRI)
{
  char bs1[32];
  char bs2[32];
  char bs3[32];
  int test;

  std::ifstream datf;
  datf.open(filename);
  datf.read((char *) nroe, sizeof(int));
  datf.read((char *) nroa, sizeof(int));
  datf.read((char *) nroao, sizeof(int));
  datf.read((char *) naux_1, sizeof(int));
  datf.read((char *) naux_2, sizeof(int));
  datf.read((char *) nrofint,   sizeof(long long int));
  datf.read((char *) nrofaux,   sizeof(long long int));
  datf.read((char *) nrofaux2,  sizeof(long long int));


  *coord   = new double[3*(*nroa)];
  *charges = new double[(*nroa)];
  *zeff    = new double[(*nroa)];
  *mass    = new double[(*nroa)];

  datf.read((char *) *coord, sizeof(double)*3*(*nroa));
  datf.read((char *) *charges, sizeof(double)*(*nroa));
  datf.read((char *) *zeff, sizeof(double)*(*nroa));
  datf.read((char *) *mass,    sizeof(double)*(*nroa));
  datf.read((char *) bs1,32);
  datf.read((char *) bs2,32);
  datf.read((char *) bs3,32);
  datf.close();

  std::stringstream tr1,tr2,tr3;
  tr1 <<bs1;
  tr1>> (*basisNameOB);
  tr2 <<bs2;
  tr2 >> (*basisNameJK);
  tr3 <<bs3;
  tr3 >> (*basisNameRI);


}


void read_oei(std::string filename,int nroao,double* Hmat,double* Tmat,double* Smat,
              double* Vmat)
{
  std::ifstream datf;
  datf.open(filename);
  datf.read((char * ) Hmat, sizeof(double)*nroao*nroao);
	datf.read((char * ) Tmat, sizeof(double)*nroao*nroao);
	datf.read((char * ) Smat, sizeof(double)*nroao*nroao);
	datf.read((char * ) Vmat, sizeof(double)*nroao*nroao);
  datf.close();
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      REM COM                                                  */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int rem_com(char* filename, char* streamstring, int string_length){
  const char com_B = '#';
  const char com_E = '\n';

  int pos = 0;
  char cc;

  std::ifstream inf(filename);

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
    std::cerr << "Buffer size exceeded !\n"; exit(0);
  }
  return(strlen(streamstring));
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*              read_inpuy (char* sysfile, int* nroao, int* nroa,             */
/*                              long long int* nrofint)                          */
/*                                                                               */
/* Read the system sizes form sysfile.                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void read_input(std::ifstream* inputfile)
{
  int nroe;
  std::string tmpline;
  std::stringstream ss_input;
  if (inputfile->is_open())
  {
    while(std::getline((*inputfile),tmpline))
    {
      //remove comments
      if (tmpline.find_first_of("#")!=std::string::npos)
        tmpline.erase(tmpline.find_first_of("#"),tmpline.length());
      if (tmpline.find_first_of(" ")!=std::string::npos)
        tmpline.erase(tmpline.find_first_of(" "),tmpline.length());
      if (tmpline.compare("\n")!=0)
        ss_input <<tmpline<<" ";
    }
  }
  inputfile->close();
  //parse input variables here
  //ss_input >> (*sysfile) >> (nroe) >> (*llim) >> (*ulim) >> (*wavefile);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                get_sys_size(char* sysfile, int* nroao, int* nroa,             */
/*                              long long int* nrofint)                          */
/*                                                                               */
/* Read the system sizes form sysfile.                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void get_sys_size(std::string sysfile,int* nroe, int* nroao, int* nroa, long long int* nrofint,char* basisName){

  std::ifstream inf(sysfile.c_str());

  inf.read((char *) nroe,    sizeof(int));
  inf.read((char *) nroao,   sizeof(int));
  inf.read((char *) nroa,    sizeof(int));
  inf.read((char *) nrofint, sizeof(long long int));
  inf.read((char *) basisName, 32);

  inf.close();
}


void read_tei(std::string filename,long long int nrofint, long long int* sortcount,double* intval,
              unsigned short* intnums)
{
  std::ifstream datf;
  datf.open(filename);
  datf.read((char *) sortcount, sizeof(long long int)*4);
	datf.read((char *) intval,    sizeof(double)*nrofint);
	datf.read((char *) intnums,   sizeof(unsigned short)*nrofint*4);
	datf.close();

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
  std::ifstream datf(sysfile.c_str());

  int nroao, nroa, nroe;
  long long int nrofint;
  char basisName[32];
  //SYSTEM DATA
  datf.read((char *) &nroe,  sizeof(int));
  datf.read((char *) &nroao , sizeof(int));
  datf.read((char *) &nroa  , sizeof(int));
  datf.read((char *) &nrofint,  sizeof(long long int));
  datf.read((char *) basisName, 32);


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
/*               write_wav_HF(........)                                          */
/*                                                                               */
/*                                                                               */
/* write  HF wave_function                                                       */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void write_wav_HF(std::string wavfile, int nroao, double* MOens, double* MOs){
  std::ofstream outf(wavfile.c_str());
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
  std::ifstream inf(wavfile.c_str());
  int real_nroao;

  inf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao != nroao){
    std::cerr << "Wrong HF wavefunction size in read_wav_HF!\n"; exit(3);
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


void   output_matrix(double* mat,int nr_of_col, int nrop, std::ofstream *outf){
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
