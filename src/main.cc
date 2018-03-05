#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "MATH/ops_mat.h"
#include "UTIL/util.h"
#include "IO/ops_io.h"
#include "SCF/ops_rhf.h"
#include "SCF/scf.h"
#include "POSTHF/motrans.h"
#include "POSTHF/loc_mp2.h"
#include "POSTHF/can_mp2.h"
#include "POSTHF/ri.h"
#include "POSTHF/ccsd.h"
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

    systeminfo *sysinfo = new systeminfo;

    std::map<std::string,std::pair<double,double> > timer;

    timer["total"] = std::make_pair(0.0,0.0);
    timer["total"].first = omp_get_wtime();


    if(argc != 2) {
        std::cerr << "Need PREFIX\n";
        exit(1);
    }

    sysinfo->prefix = argv[1];
    read_system(sysinfo);

    sysinfo->ion_rep =  calc_ion_rep(sysinfo);
    print_header(sysinfo);


    libint2::initialize();
    libint2::BasisSet obs(sysinfo->basisNameOB,sysinfo->atoms);
    libint2::BasisSet dfbs(sysinfo->basisNameRI,sysinfo->atoms);
    obs.set_pure(true);dfbs.set_pure(true);

    /*
     *
     * PART2 : One Electron & Electron Matrices and Integrals
     *
     */

    OEints* onemats = new OEints;
    TEints* twomats = new TEints;

    allocate_onemats(sysinfo,onemats);
    // Always read the One-Electron Matrices
    // Maybe on a day in a far far future, we calculate the OEI

    read_oei(sysinfo,onemats);


    //Two possibilities:
    //1) ERI provided -> read them
    //2) NO ERI provided -> Transform Hmat to Libint -> Calculate -> ERI

    if (sysinfo->nrofint > 0) {
        std::cout << "Two electron integrals provided, reading in .." << '\n';
        twomats->intval        = new double[sysinfo->nrofint];                   //two electron integrals
        twomats->intnums       = new unsigned short[sysinfo->nrofint*4];         //two electron indices

        //ints_start = omp_get_wtime();
        read_tei(sysinfo,twomats);
        //ints_end = omp_get_wtime();
    }else{
        std::cout << "No Tei's provided ... using libint to calculate them" << '\n';

        calculate_libint_oei(sysinfo,obs,onemats);
        //ints_start = omp_get_wtime();
        calculate_libint_tei(sysinfo,obs,twomats);
        //ints_end = omp_get_wtime();

    }
    calc_S12(sysinfo,onemats);

    /*
     *
     * Part 3: Read Guess & Perform SCF
     *
     */


    read_wav_HF(sysinfo,onemats);
    check_initial_guess(sysinfo,onemats);

    sysinfo->scfiter = 100;
    timer["SCF"] = std::make_pair(0.0,0.0);
    timer["SCF"].first =omp_get_wtime();
    run_scf(sysinfo,onemats,twomats);
    timer["SCF"].second=omp_get_wtime();


    /*
     * PART 4:
     * Post Hartree-Fock Part:
     * BPQ of /c basis
     * transform to Bia
     * perform mp2
     *
     */

    pHF* postHF = new pHF;

    build_FMo(sysinfo,onemats,postHF);
    {

        double trafo_start = omp_get_wtime();
        MOtrans(sysinfo,onemats,twomats,postHF);

        double trafo_end = omp_get_wtime();
        double mp2_start = omp_get_wtime();
        run_canonical_mp2(sysinfo,postHF);
        double mp2_end   = omp_get_wtime();
    }



    //RI-MP2 //trafo

    postHF->BPQ = new double[sysinfo->naux_2*sysinfo->nroao*sysinfo->nroao];
    calculate_ri(obs,dfbs,postHF->BPQ);
    double ritrafo_start = 0.0;
    double ritrafo_end   = 0.0;
    double rimp2_start   = 0.0;
    double rimp2_end     = 0.0;

    {
        ritrafo_start = omp_get_wtime();
        postHF->Bia = new double[(sysinfo->naux_2)*(sysinfo->nroe/2)*(sysinfo->nroao-sysinfo->nroe/2)];
        //read_transform_ri(prefix,nroe,nroao,naux_2,MOs,Bia);
        transform_ri(sysinfo,onemats,postHF);
        ritrafo_end = omp_get_wtime();
        rimp2_start = omp_get_wtime();
        run_canonical_mp2_ri(sysinfo,postHF);
        rimp2_end = omp_get_wtime();
     } //end RI-MP2


    ccsd_ur(sysinfo,onemats,postHF);


    /* What the fuck!! do we need ?
     * FOCK Matrix in spin-orbital basis
     * t1 - amplides
     * t2 - amplitudes
     * eri tensor in spin orbital basis
     */

/*
    std::cout<<"nel:"<<nroe<<std::endl;
    std::cout<<"nocc:"<<nroe<<std::endl;
    std::cout<<"nvir:"<<(2*nroao - nroe)<<std::endl;
    std::cout<<"nroao:"<<nroao<<std::endl;
    std::cout<<"nob:"<<2*nroao<<std::endl;
    std::cout<<"nalpha:"<<nroao<<std::endl;
    std::cout<<"nbeta:"<<nroao<<std::endl;




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
