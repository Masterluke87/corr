#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <libint2.hpp>

struct pHF{
    double* FMo;
    double* prec_ints;
    double* BPQ;
    double* Bia;
};


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

#endif

