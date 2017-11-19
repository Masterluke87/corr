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
void get_sys_size(std::string sysfile,int* nroe, int* nroao, int* nroa, long long int* nrofint);
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
