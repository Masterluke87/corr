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



void print_header();
int  rem_com(char* filename, char* streamstring, int string_length);
void get_sys_size(std::string sysfile,int* nroe, int* nroao, int* nroa, long long int* nrofint,char* basisName);
void read_sys(std::string sysfile, double* coord, double* charges, double* mass,
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy,
	      double *Dz, long long int* sortcount, double* intval,
	      unsigned short* intnums);
void read_sys_1el(char* sysfile, double* coord, double* charges, double* mass,
		  double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy,
		  double *Dz);
void write_wav_HF(std::string wavfile, int nroao, double* MOens, double* MOs);
void read_wav_HF(std::string wavfile, int nroao, double* MOens, double* MOs);
void output_matrix(double* mat,int nr_of_col, int nrop, std::ofstream *outf);

void read_system(std::string filename, int* nroe,int* nroa,
                  int* nroao, int* naux_1, int* naux_2,
                  long long int* nrofint,long long int* nrofaux,long long int* nrofaux2,
                  double** coord,double** charges,double** zeff,double** masses,
                  std::string* basisNameOB, std::string* basisNameJK, std::string* basisNameRI);

void read_oei(std::string filename,int nroao,double* Hmat,double* Tmat,double* Smat,double* Vmat);
void read_tei(std::string filename,long long int nrofint,long long int* sortcount,double* intval,
              unsigned short* intnums);
