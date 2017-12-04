#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "IO/ops_io.h"
#include "SCF/ops_rhf.h"
#include "SCF/scf.h"
#include "POSTHF/motrans.h"
#include "POSTHF/loc_mp2.h"
#include "POSTHF/can_mp2.h"
#include <iomanip>
#include <omp.h>
#include <libint2.hpp>

#define PWIDTH_L 12
#define PWIDTH_R 16



int main(int argc, char const *argv[]) {
    int worldsize = 1;
    int rank = 0;

    int nroe;                  //Nr of electrons  (if negative, read in center of mass)
    int nroao;
    int nroa;
    int naux_1;
    int naux_2;

    long long int nrofint;     //Nr of two electron Integrals
    long long int nrofaux;
    long long int nrofaux2;

    double*        coord;       //atomic coordinats               3*nroa
    double*        charges;     //atomic charges                    nroa
    double*              zeff;              //in case ecps are used
    double*        mass;        //atomic masses                     nroa

    std::string basisNameOB;
    std::string basisNameJK;
    std::string basisNameRI;

    double total_start = omp_get_wtime();
    print_header();

    if(argc != 2) {
        std::cerr << "Need PREFIX\n";
        exit(1);
    }

    std::string prefix = argv[1];

    read_system(prefix+".sys",&nroe,&nroa,&nroao,&naux_1,
                &naux_2,&nrofint,&nrofaux,&nrofaux2,&coord,
                &charges,&zeff,&mass,&basisNameOB,&basisNameJK,&basisNameRI);

    double ion_rep =  calc_ion_rep( nroa, coord, zeff);

    {  //output

        std::cout << "\nSYSTEMDATA" << '\n';
        std::cout << "----------" << "\n\n";
        std::cout<<std::setw(-2)<<"Z" <<std::setw(10)<<"zeff"<<std::setw(10)<<"x" <<std::setw(10) << "y" <<std::setw(10)<< "z" <<'\n';
        for (size_t i = 0; i < nroa; i++) {
            std::cout<<std::setw(-2)<<charges[i]<<std::setw(10)<<zeff[i]<<std::setw(10)<< coord[i*3+0]<<std::setw(10) << coord[i*3+1]<<std::setw(10)<< coord[i*3+2]<<'\n';
        }

        std::cout << "nroe          :" << nroe<<'\n';
        std::cout << "nroa          :" << nroa<<'\n';
        std::cout << "nroao         :" << nroao<<'\n';
        std::cout << "naux-JK       :" << naux_1<<'\n';
        std::cout << "naux-RI       :" << naux_2<<'\n';
        std::cout << "nrofint       :" << nrofint<<'\n';
        std::cout << "nrofaux(JK)   :" << nrofaux<<'\n';
        std::cout << "nrofaux(RI)   :" << nrofaux2<< '\n';
        std::cout << "basisName(OB) :" << basisNameOB<<'\n';
        std::cout << "basisName(JK) :" << basisNameJK<<'\n';
        std::cout << "basisName(RI) :" << basisNameRI<<'\n';
    }

    double scf_start;
    double ints_start;

    double scf_end;
    double ints_end;



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

    double*        Tmat_libint;
    double*        Hmat_libint;
    double*              Smat_libint;
    double*              Vmat_libint;

    double*        Tmat_trans;
    double*        Hmat_trans;
    double*        Smat_trans;
    double*        Vmat_trans;

    dumd = new double[17*nroao*nroao];
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
    Tmat_libint = &(dumd[inc]); inc+=nroao*nroao;
    Hmat_libint = &(dumd[inc]); inc+=nroao*nroao;

    Smat_libint = &(dumd[inc]); inc+=nroao*nroao;
    Vmat_libint = &(dumd[inc]); inc+=nroao*nroao;
    Tmat_trans = &(dumd[inc]); inc+=nroao*nroao;
    Hmat_trans = &(dumd[inc]); inc+=nroao*nroao;
    Smat_trans = &(dumd[inc]); inc+=nroao*nroao;

    Vmat_trans = &(dumd[inc]); inc+=nroao*nroao;

    read_oei(prefix+".oei",nroao,Hmat,Tmat, Smat,Vmat);



    double* intval;
    unsigned short* intnums;
    long long int sortcount[4];

    if (nrofint > 0) {
        std::cout << "Two electron integrals provided, reading in .." << '\n';
        intval        = new double[nrofint];                   //two electron integrals
        intnums       = new unsigned short[nrofint*4];         //two electron indices

        ints_start = omp_get_wtime();
        read_tei(prefix+".tei",nrofint,sortcount,intval,intnums);
        ints_end = omp_get_wtime();

        calc_S12(nroao, Smat, Som12);

        read_wav_HF(prefix+".ahfw",nroao,MOens,MOs);
        { //check of gd orbitals are provided;
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
        scf_start = omp_get_wtime();
        run_scf(nroao,nroe,MOs,Pmat,Hmat,Fmat,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);
        scf_end = omp_get_wtime();
    }else{
        std::cout << "No Tei's provided ... using libint to calculate them" << '\n';
        libint2::initialize();

        std::cout << "Starting libint2 - part" << '\n';
        std::vector<libint2::Atom> atoms(nroa);
        for (size_t i = 0; i < nroa; i++) {
            atoms[i].atomic_number = charges[i];
            atoms[i].x = coord[i*3+0];
            atoms[i].y = coord[i*3+1];
            atoms[i].z = coord[i*3+2];
        }
        libint2::BasisSet obs(basisNameOB,atoms);
        obs.set_pure(true);
        calculate_libint_oei(atoms,obs,zeff,
                             Hmat,        Tmat, Smat, Vmat,
                             Hmat_libint, Tmat_libint, Smat_libint, Vmat_libint,
                             Hmat_trans, Tmat_trans, Smat_trans, Vmat_trans);
        Smat = Smat_trans;
        Hmat = Hmat_trans;

        ints_start = omp_get_wtime();
        calculate_libint_tei(atoms,obs,nrofint,&intval,&intnums,sortcount);
        ints_end = omp_get_wtime();
        std::cout<<"Integral evaluation: "<<std::setw( 4 )<<ints_start - ints_end <<" [s]"<<"\n\n";
        calc_S12(nroao, Smat, Som12);
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
        scf_start = omp_get_wtime();
        run_scf(nroao,nroe,MOs,Pmat,Hmat,Fmat,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);
        scf_end = omp_get_wtime();
        libint2::finalize();
    }
/*
    //RI-with libint
    {
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
     libint2::Engine eri2_engine(libint2::Operator::coulomb, dfbs.max_nprim(), dfbs.max_l());
     eri2_engine.set_braket(libint2::BraKet::xs_xs);

     auto shell2bf = dfbs.shell2bf();
     const auto& buf_vec_eri = eri2_engine.results();

     double * J = new double[dfbs.nbf()*dfbs.nbf()];
     int naux = dfbs.nbf();

     for(auto s2=0; s2!=dfbs.size(); ++s2)
         for(auto s1=s2; s1!=dfbs.size(); ++s1){
             eri2_engine.compute(dfbs[s1], dfbs[s2]);
             auto ints_shellset_eri = buf_vec_eri[0];
             if (ints_shellset_eri!=NULL){

                 auto bf1 = shell2bf[s1]; // first basis function in first shell
                 auto n1 = dfbs[s1].size(); // number of basis functions in first shell
                 auto bf2 = shell2bf[s2]; // first basis function in second shell
                 auto n2 = dfbs[s2].size(); // number of basis functions in second shell

                 for(auto f1=0; f1!=n1; ++f1)
                     for(auto f2=0; f2!=n2; ++f2)
                         J[(bf1+f1)*naux + (bf2+f2)] = ints_shellset_eri[f1*n2+f2];
            }
         }
    }
    */
    libint2::finalize();





    long long int prec_mem = (long long int) nroe/2*(long long int) nroao* (long long int) nroao* (long long int) nroao;
    std::cout << "Need " << prec_mem*sizeof(double) << " bytes (" << prec_mem*sizeof(double)/1024/1024 << " MB) for precalculation\n";
    std::cout << prec_mem << " doubles on "<< worldsize<< " cpus \n";
    if (worldsize > 1) {
        std::cout << "Allocating " << prec_mem/worldsize << " on "<< (worldsize-1)<< " processors and "<< (prec_mem/worldsize + prec_mem%worldsize)<< " on one cpu";
    }
    double* prec_ints;
    if (rank==worldsize-1) {
        prec_ints =  new double[prec_mem/worldsize + prec_mem%worldsize];
    }else{
        prec_ints =  new double[prec_mem/worldsize];
    }




    //4-index trans
    double trafo_start = omp_get_wtime();
    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;
    int i,j,k,l;
    long long int prec_count = 0;
    for(i = 0; i < nroe/2; i++) {
        for(j = 0; j < nroao; j++) {
            for(k = 0; k < nroao; k++) {
                for(l = 0; l <  nroao; l++) {
                    prec_ints[prec_count%(prec_mem/worldsize)] = mo2int_op(i, j, k, l,
                                                                           nroao, MOs, nrofint, sortcount, intval, intnums,&std::cout);
                    prec_count++;
                }
            }
        }
        std::cout << i << "\t"<<std::flush;
        if((i+1)%10==0) std::cout << "\n"<<std::flush;
    }

    double trafo_end = omp_get_wtime();
    //
    //calculate the Hartree-Fock energy:
    //
    double E_hf = 0.0;
    double oneE = 0.0;
    double twoE = 0.0;
    int mu = 0;
    int nu = 0;
    double* HMo = new double[nroao*nroao];
    double* GMo = new double[nroao*nroao];
    double* FMo = new double[nroao*nroao];
    for (size_t i = 0; i < nroao*nroao; i++) {
        HMo[i] = 0.0;
        GMo[i] = 0.0;
        FMo[i] = 0.0;
    }

    for (int i = 0; i < nroao; i++) {
        for (int j = 0; j < nroao; j++) {
            for (int mu = 0; mu < nroao; mu++) {
                for (int nu = 0; nu < nroao; nu++) {
                    HMo[i*nroao + j] += MOs[i*nroao + mu] * Hmat[mu*nroao + nu] * MOs[j*nroao + nu];

                }
            }
        }
    }

    for (int i = 0; i < nroao; i++) {
        for (int j = 0; j < nroao; j++) {
            for (int k = 0; k < nroe/2; k++) {
                GMo[i*nroao + j] += (2*prec_ints[k*istep + k*jstep + i*kstep + j]
                                     -prec_ints[k*istep + i*jstep + k*kstep + j]);
            }
        }
    }

    for (int i = 0; i < nroao; i++) {
        for (int j = 0; j < nroao; j++) {
            FMo[i*nroao + j] = HMo[i*nroao + j] + GMo[i*nroao + j];
        }
    }

    //    run_non_canonical_mp2(nroao,nroe,prec_ints,FMo);

    for (size_t i = 0; i < nroe/2; i++) {
        oneE += 2 * HMo[i*nroao + i];
    }
    oneE += ion_rep;
    std::cout<<"\n\n\n"<<std::setw( PWIDTH_L )<< "oneE:" <<std::setw( 16 )<<oneE<<'\n';

    for (size_t i = 0; i < nroe/2; i++) {
        for (size_t j = 0; j < nroe/2; j++) {
            twoE+= (2 *prec_ints[i*istep + i*jstep + j*kstep + j]
                    - prec_ints[i*istep + j*jstep + i*kstep + j]);

        }
    }
    std::cout << std::setw( PWIDTH_L ) << "twoE:" <<std::setw( PWIDTH_R )<<twoE<< "\n";
    std::cout << std::setw( PWIDTH_L ) << "HF:"   <<std::setw( PWIDTH_R )<<oneE+twoE<<"\n";
    ///END Hartree-Fock


    //RI-MP2
    double ritrafo_start = 0.0;
    double ritrafo_end   = 0.0;
    {
        ritrafo_start = omp_get_wtime();
        int nocc = nroe/2;
        int nvir = nroao-nroe/2;
        double* BPQ = new double[naux_2*nroao*nroao];
        double* BiQ = new double[naux_2*nroao*nroao];
        double* Bia = new double[naux_2*nroao*nroao];

        for(int i=0;i<naux_2*nroao*nroao;i++){
            BiQ[i] =0.0;
            Bia[i] =0.0;
        }

        std::ifstream datf;
        datf.open(prefix+".rimp2");
        datf.read((char*) BPQ , naux_2*nroao*nroao*sizeof(double));
        datf.close();

        /*
        int a,b,c,d;
        double aoint = 0.0;
        double riint = 0.0;
        double diff=0.0;
        double maxdiff=0.0;
        double cumdiff=0.0;

        for (int j=0;j<nrofint;j++)
        {
        a =  intnums[j*4+0];
        b =  intnums[j*4+1];
        c =  intnums[j*4+2];
        d =  intnums[j*4+3];
        aoint = intval[j];
        riint = 0.0;
        for (int i=0;i<naux_2;i++){
            riint += BPQ[i*nroao*nroao+a*nroao+b]*BPQ[i*nroao*nroao+c*nroao+d];
        }
        diff = std::fabs(aoint-riint);
        cumdiff+=diff;
        if (diff>maxdiff)
            maxdiff=diff;
        }
        std::cout<<"Difference(cum): "<<cumdiff;
        std::cout<<"Difference(max): "<<maxdiff;
        */
        for (int Q=0; Q<naux_2;Q++)
            for (int i=0;i<nroao;i++)
                for (int mu=0; mu<nroao; mu++)
                    for (int nu=0; nu<nroao;nu++)
                        BiQ[Q*nroao*nroao + i*nroao + nu] +=  MOs[i*nroao+mu]*BPQ[Q*nroao*nroao + mu*nroao + nu];

        for (int Q=0; Q<naux_2;Q++)
           for (int a=0;a<nroao;a++)
               for (int i =0; i<nroao;i++)
                    for (int nu=0; nu<nroao;nu++)
                        Bia[Q*nroao*nroao+i*nroao+a]+= MOs[a*nroao + nu]*BiQ[Q*nroao*nroao + i*nroao + nu ];

        ritrafo_end = omp_get_wtime();

        double iajb = 0.0;
        double ibja = 0.0;
        double EMP2,EMP2_SS,EMP2_OS;
        EMP2=EMP2_SS=EMP2_OS=0.0;
        for (int i = 0; i < nroe/2; i++) {
          for (int j = 0; j < nroe/2; j++) {
            for (int a = nroe/2; a < nroao; a++) {
              for (int b = nroe/2;  b < nroao; b++) {
                  iajb = ibja = 0.0;
                  for (int Q=0;Q<naux_2;Q++){
                      ibja += Bia[Q*nroao*nroao+i*nroao+b]*Bia[Q*nroao*nroao+j*nroao+a];
                      iajb += Bia[Q*nroao*nroao+i*nroao+a]*Bia[Q*nroao*nroao+j*nroao+b];
                  }

                EMP2_SS += - (iajb*(iajb))/ (FMo[a*nroao +a]+FMo[b*nroao+b]-FMo[i*nroao+i]-FMo[j*nroao+j]);
                EMP2_OS += - ((iajb - ibja)*iajb)/(FMo[a*nroao+a]+FMo[b*nroao+b]-FMo[i*nroao +i]-FMo[j*nroao + j]);
              }
            }
          }
        }

        EMP2 = EMP2_OS + EMP2_SS;
        std::cout<<"\nRI-Results:";
        std::cout <<std::setw( PWIDTH_L ) << "EMP2_SS:" <<std::setw( PWIDTH_R )<<EMP2_SS<< '\n';
        std::cout <<std::setw( PWIDTH_L ) << "EMP2_OS:" <<std::setw( PWIDTH_R )<<EMP2_OS<< '\n';
        std::cout <<std::setw( PWIDTH_L ) << "EMP2:"    <<std::setw( PWIDTH_R )<<EMP2<< '\n';


        /* DEBUG
        long long int kstep = nroao;
        long long int jstep = nroao*kstep;
        long long int istep = nroao*jstep;

        int i,j,k,l;
        i=j=k=l=0;
        for (int i=0; i<nocc;i++)
            for (int j=0;j<nroao;j++)
                for (int k=0;k<nroao;k++)
                    for (int l=0;l<nroao;l++){
                riint =0.0;
            for (int Q=0;Q<naux_2;Q++)
            {
                riint += Bia[Q*nroao*nroao+i*nroao+j] *Bia[Q*nroao*nroao+k*nroao+l];
            }
           // std::cout<< prec_ints[i*istep+j*jstep+k*kstep+l] << " "<<riint<<"\n";
        } */
        delete[] BPQ;
        delete[] BiQ;
        delete[] Bia;

    }


    double mp2_start = omp_get_wtime();
    run_canonical_mp2(nroe,nroao,prec_ints,FMo);
    double mp2_end   = omp_get_wtime();

    std::cout <<"\n\nTIMINGS:\n";
    std::cout << std::setw( PWIDTH_L ) << "TEI [s]:" <<std::setw( PWIDTH_R )<<ints_end - ints_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "SCF [s]:" <<std::setw( PWIDTH_R )<<scf_end-scf_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "4-index [s]:" <<std::setw( PWIDTH_R )<<trafo_end-trafo_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "ri-4ind [s]:" <<std::setw( PWIDTH_R )<<ritrafo_end-ritrafo_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "MP2 [s]:" <<std::setw( PWIDTH_R )<<mp2_end-mp2_start<< "s \n";

    double total_end = omp_get_wtime();
    std::cout << "---------------------------" << '\n';
    std::cout<< std::setw( PWIDTH_L ) << "Total [s]:"<<std::setw( PWIDTH_R )<<total_end-total_start<<"s \n";
    std::cout << "\n\n--END--" << '\n';

    return 0;
}
