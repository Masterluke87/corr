#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "MATH/ops_mat.h"
#include "IO/ops_io.h"
#include "SCF/ops_rhf.h"
#include "SCF/scf.h"
#include "POSTHF/motrans.h"
#include "POSTHF/loc_mp2.h"
#include "POSTHF/can_mp2.h"
#include "POSTHF/ri.h"
#include <iomanip>
#include <omp.h>
#include <libint2.hpp>
#include <cblas.h>
#include <lapacke.h>


#define PWIDTH_L 12
#define PWIDTH_R 16



int main(int argc, char const *argv[]) {
    int worldsize = 1;
    int rank = 0;

    /*
     *
      Part 1 : Initialization
     *
     */

    // System Information:

    int nroe;                  //Nr of electrons  (if negative, read in center of mass)
    int nroa;                  //nr of atoms

    int nroao;                 //nr of basis functions
    int naux_1;                //nr of /JK
    int naux_2;                //nr of /C

    std::string basisNameOB;
    std::string basisNameJK;
    std::string basisNameRI;


    double*        coord;       //atomic coordinats               3*nroa
    double*        charges;     //atomic charges                    nroa
    double*        zeff;        //in case ecps are used             nroa
    double*        mass;        //atomic masses                     nroa

    long long int nrofint;     //Nr of two electron Integrals
    long long int nrofaux;
    long long int nrofaux2;

    std::map<std::string,std::pair<double,double> > timer;

    timer["total"] = std::make_pair(0.0,0.0);
    timer["total"].first = omp_get_wtime();


    if(argc != 2) {
        std::cerr << "Need PREFIX\n";
        exit(1);
    }

    std::string prefix = argv[1];
    read_system(prefix+".sys",&nroe,&nroa,&nroao,&naux_1,
                &naux_2,&nrofint,&nrofaux,&nrofaux2,&coord,
                &charges,&zeff,&mass,&basisNameOB,&basisNameJK,&basisNameRI);
    double ion_rep =  calc_ion_rep( nroa, coord, zeff);
    print_header(&nroe,&nroa,&nroao,&naux_1,&naux_2,&nrofint,&nrofaux,
                 &nrofaux2,coord,charges,zeff,
                 &basisNameOB,&basisNameJK,&basisNameRI);


    libint2::initialize();
    std::vector<libint2::Atom> atoms(nroa);
    for (size_t i = 0; i < nroa; i++) {
        atoms[i].atomic_number = charges[i];
        atoms[i].x = coord[i*3+0];
        atoms[i].y = coord[i*3+1];
        atoms[i].z = coord[i*3+2];
    }
    libint2::BasisSet obs(basisNameOB,atoms);
    libint2::BasisSet dfbs(basisNameRI,atoms);
    obs.set_pure(true);
    dfbs.set_pure(true);

    /*
     *
     * PART2 : One Electron & Electron Matrices and Integrals
     *
     */



    double* dumd;
    //one electron mat&vecs
    double*        Hmat;        //one electron Hamiltionian         nroao*nroao
    double*        Tmat;        //Kinetic energy operator           nroao*nroao
    double*        Smat;        //Overlap matrix S                  nroao*nroao
    double*        Som12;       //S^-1/2                            nroao*nroao
    double*        Vmat;        //Nuclear repuslsion

    double*        MOens;       //MO Energies
    double*        MOs;         //MO coeffs
    double*        Pmat;      //density matrix                    nroao*nroao
    double*        Fmat;

    dumd = new double[9*nroao*nroao];
    int inc = 0;

    Hmat  = &(dumd[inc]); inc+=nroao*nroao;
    Tmat  = &(dumd[inc]); inc+=nroao*nroao;
    Smat  = &(dumd[inc]); inc+=nroao*nroao;
    Som12 = &(dumd[inc]); inc+=nroao*nroao;
    Vmat  = &(dumd[inc]); inc+=nroao*nroao;

    MOens = &(dumd[inc]); inc+=nroao*nroao;
    MOs   = &(dumd[inc]); inc+=nroao*nroao;
    Pmat  = &(dumd[inc]); inc+=nroao*nroao;
    Fmat  = &(dumd[inc]); inc+=nroao*nroao;


    // Always read the One-Electron Matrices
    // Maybe on a day in a far far future, we calculate the OEI

    read_oei(prefix+".oei",nroao,Hmat,Tmat, Smat,Vmat);
    double* intval;
    unsigned short* intnums;
    long long int sortcount[4];

    //Two possibilities:
    //1) ERI provided -> read them
    //2) NO ERI provided -> Transform Hmat to Libint -> Calculate -> ERI

    if (nrofint > 0) {
        std::cout << "Two electron integrals provided, reading in .." << '\n';
        intval        = new double[nrofint];                   //two electron integrals
        intnums       = new unsigned short[nrofint*4];         //two electron indices

        //ints_start = omp_get_wtime();
        read_tei(prefix+".tei",nrofint,sortcount,intval,intnums);
        //ints_end = omp_get_wtime();
    }else{
        std::cout << "No Tei's provided ... using libint to calculate them" << '\n';

        calculate_libint_oei(atoms,obs,zeff,Hmat,Tmat, Smat, Vmat);
        //ints_start = omp_get_wtime();
        calculate_libint_tei(atoms,obs,nrofint,&intval,&intnums,sortcount);
        //ints_end = omp_get_wtime();

    }
    calc_S12(nroao, Smat, Som12);

    /*
     *
     * Part 3: Read Guess & Perform SCF
     *
     */
    read_wav_HF(prefix+".ahfw",nroao,MOens,MOs);
    {             //check if gd orbitals are provided;
        double modiag=0.0;
        for (size_t i = 0; i < nroao; i++) {
            modiag += MOs[i*nroao+i];
        }
        if (std::fabs(modiag) < 0.1) {
            std::cout << "No converged MOs provided... doing core guess" << '\n';
            for (int i =0; i<nroao*nroao; i++)
                Fmat[i] = Hmat[i];
            double *tmpmat = new double[nroao*nroao];
            diag_Fmat(nroao, Fmat,MOs,MOens,Som12, tmpmat);
            delete[] tmpmat;
        }
    }
    timer["SCF"] = std::make_pair(0.0,0.0);
    timer["SCF"].first =omp_get_wtime();
    run_scf(nroao,nroe,MOs,Pmat,Hmat,Smat,Fmat,MOens,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);
    timer["SCF"].second=omp_get_wtime();


    /*
     * PART 4:
     * Post Hartree-Fock Part:
     * BPQ of /c basis
     * transform to Bia
     * perform mp2
     *
     */

    double* FMo  = new double[nroao*nroao];
    build_FMo(nroao,Fmat,FMo,MOs);


    {
        double* prec_ints;
        double trafo_start = omp_get_wtime();
        MOtrans(MOs,nroao,nroe,nrofint,sortcount,intval,intnums,&prec_ints);
        double trafo_end = omp_get_wtime();

        double mp2_start = omp_get_wtime();
        run_canonical_mp2(nroe,nroao,prec_ints,FMo);
        double mp2_end   = omp_get_wtime();
    }



    //RI-MP2 //trafo

    double* BPQ = new double[naux_2*nroao*nroao];
    calculate_ri(obs,dfbs,BPQ);
    double ritrafo_start = 0.0;
    double ritrafo_end   = 0.0;
    double rimp2_start   = 0.0;
    double rimp2_end     = 0.0;
    {
        ritrafo_start = omp_get_wtime();
        double* Bia = new double[naux_2*nroe/2*(nroao-nroe/2)];
        //read_transform_ri(prefix,nroe,nroao,naux_2,MOs,Bia);
        transform_ri(nroe,nroao,naux_2,MOs,BPQ,Bia);
        ritrafo_end = omp_get_wtime();
        rimp2_start = omp_get_wtime();
        run_canonical_mp2_ri(nroe,nroao,naux_2,Bia,FMo);
        rimp2_end = omp_get_wtime();
        delete[] Bia;
    } //end RI-MP2




    /* Part 5
     * Print timings
     */



    std::cout <<"\n\nTIMINGS:\n";
    std::cout <<"--------\n";
   // std::cout << std::setw( PWIDTH_L ) << "TEI [s]:" <<std::setw( PWIDTH_R )<<ints_end - ints_start<< " s \n";
   // std::cout << std::setw( PWIDTH_L ) << "SCF [s]:" <<std::setw( PWIDTH_R )<<scf_end-scf_start<< " s \n";
   // std::cout << std::setw( PWIDTH_L ) << "4-index [s]:" <<std::setw( PWIDTH_R )<<trafo_end-trafo_start<< " s \n";
   // std::cout << std::setw( PWIDTH_L ) << "ri-4ind [s]:" <<std::setw( PWIDTH_R )<<ritrafo_end-ritrafo_start<< " s \n";
   // std::cout << std::setw( PWIDTH_L ) << "RI-MP2 [s]:" <<std::setw( PWIDTH_R )<<rimp2_end-rimp2_start<< " s \n";
   // std::cout << std::setw( PWIDTH_L ) << "MP2 [s]:" <<std::setw( PWIDTH_R )<<mp2_end-mp2_start<< " s \n";

    timer["total"].second = omp_get_wtime();
    std::cout << "------------------------------" << '\n';
    std::cout<< std::setw( PWIDTH_L ) << "Total [s]:"<<std::setw( PWIDTH_R )<<timer["total"].second - timer["total"].first<<" s \n";
    std::cout << "\n\n--END--" << '\n';
    libint2::finalize();
    return 0;
}
