#ifndef OPS_IO_H
#define OPS_IO_H


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
#include <libint2.hpp>


struct systeminfo{

    int nroe;
    int nroao;
    int nroa;

    int naux_1;
    int naux_2;

    long long int nrofint;     //Nr of two electron Integrals
    long long int nrofaux;
    long long int nrofaux2;


    std::string basisNameOB;
    std::string basisNameJK;
    std::string basisNameRI;

    double*        coord;       //atomic coordinats               3*nroa
    double*        charges;     //atomic charges                    nroa
    double*        zeff;        //in case ecps are used             nroa
    double*        mass;        //atomic masses                     nroa

    std::vector<libint2::Atom> atoms;

    double ion_rep;
    int scfiter;

    std::string prefix;
};

struct OEints{
    double*        Hmat;        //one electron Hamiltionian         nroao*nroao
    double*        Tmat;        //Kinetic energy operator           nroao*nroao
    double*        Smat;        //Overlap matrix S                  nroao*nroao
    double*        Som12;       //S^-1/2                            nroao*nroao
    double*        Vmat;        //Nuclear repuslsion

    double*        MOens;       //MO Energies
    double*        MOs;         //MO coeffs
    double*        Pmat;        //density matrix                    nroao*nroao
    double*        Fmat;        //Fock-Matrix

};

struct TEints{
    double* intval;
    unsigned short* intnums;
    long long int sortcount[4];
};



void print_header(systeminfo *sysinfo);

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
void read_wav_HF(systeminfo *sysinfo, OEints * onemats);
void output_matrix(double* mat,int nr_of_col, int nrop, std::ofstream *outf);

void read_system(systeminfo* sysinfo);
void read_oei(systeminfo *sysi, OEints *oemats);
void read_tei(systeminfo *sysi, TEints *t);


#endif
