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
    double*        zeff;        //in case ecps are used             nroa
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

    {     //output

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

    std::map<std::string,std::pair<double,double> > timer;



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
    double*        Smat_libint;
    double*        Vmat_libint;

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

        ints_start = omp_get_wtime();
        read_tei(prefix+".tei",nrofint,sortcount,intval,intnums);
        ints_end = omp_get_wtime();
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
        libint2::finalize();
    }
    calc_S12(nroao, Smat, Som12);
    read_wav_HF(prefix+".ahfw",nroao,MOens,MOs);
    {             //check of gd orbitals are provided;
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
    scf_start = omp_get_wtime();
    run_scf(nroao,nroe,MOs,Pmat,Hmat,Smat,Fmat,MOens,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);
    timer["SCF"].second=omp_get_wtime();
    scf_end = omp_get_wtime();



    //RI-with libint
    double* BPQ;

    {
        libint2::initialize();
        std::vector<libint2::Atom> atoms(nroa);
        for (size_t i = 0; i < nroa; i++) {
            atoms[i].atomic_number = charges[i];
            atoms[i].x = coord[i*3+0];
            atoms[i].y = coord[i*3+1];
            atoms[i].z = coord[i*3+2];
        }
        //libint2::Engine::
        libint2::BasisSet obs(basisNameOB,atoms);
        libint2::BasisSet dfbs(basisNameRI,atoms);

        libint2::Engine eri2_engine(libint2::Operator::coulomb, dfbs.max_nprim(), dfbs.max_l());
        eri2_engine.set_braket(libint2::BraKet::xs_xs);

        auto shell2bf = dfbs.shell2bf();
        const auto& buf_vec_eri = eri2_engine.results();

        if (dfbs.nbf()!=naux_2) {
            std::cout<<"DANGER!! STH IS WRONG";
            exit(0);
        }
        double * J = new double[dfbs.nbf()*dfbs.nbf()];
        memset(J,0,naux_2*naux_2*sizeof(double));


        for(auto s1=0; s1!=dfbs.size(); ++s1)
            for(auto s2=0; s2!=dfbs.size(); ++s2) {
                eri2_engine.compute(dfbs[s1], dfbs[s2]);
                auto ints_shellset_eri = buf_vec_eri[0];
                if (ints_shellset_eri!=NULL) {

                    auto bf1 = shell2bf[s1]; // first basis function in first shell
                    auto n1 = dfbs[s1].size(); // number of basis functions in first shell
                    auto bf2 = shell2bf[s2]; // first basis function in second shell
                    auto n2 = dfbs[s2].size(); // number of basis functions in second shell

                    for(auto f1=0; f1!=n1; ++f1)
                        for(auto f2=0; f2!=n2; ++f2)
                            J[(bf1+f1)*naux_2 + (bf2+f2)] = ints_shellset_eri[f1*n2+f2];
                }
            }

        double* eigval   = new double[naux_2];
        double* eigvec   = new double[naux_2*naux_2];
        double* eigvec_c = new double[naux_2*naux_2];


        std::memcpy(eigvec,J,naux_2*naux_2*sizeof(double));

        LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',naux_2,eigvec,naux_2,eigval);

        //debug
        double* tmpJ = new double[naux_2*naux_2];
        memset(tmpJ,0,naux_2*naux_2*sizeof(double));
        {
            for(int z=0;z<naux_2;z++)
                for(int y=0;y<naux_2;y++)
                    for(int x=0;x<naux_2;x++)
                        tmpJ[y*naux_2+z] += eigvec[y*naux_2 +x]* eigval[x] * eigvec[z*naux_2+x];
        }








        std::memcpy(eigvec_c,eigvec,naux_2*naux_2*sizeof(double));

        double max_eig = eigval[naux_2-1];
        double tol = 1.0E-10;
        for (int i = 0; i<naux_2; i++) {
            if(eigval[i]/max_eig < tol || eigval[i] <0)
                eigval[i] =0.0;
            else
                eigval[i] = 1 / std::sqrt(eigval[i]);
          //  cblas_dscal(naux_2,eigval[i],&eigvec[i*naux_2],1);
        }

//        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,naux_2,naux_2,
//                    naux_2,1.0,eigvec_c,naux_2,eigvec,naux_2,0.0,J,naux_2);
        {
            memset(J,0,naux_2*naux_2*sizeof(double));
            for(int z=0;z<naux_2;z++)
                for(int y=0;y<naux_2;y++)
                    for(int x=0;x<naux_2;x++)
                        J[y*naux_2+z] += eigvec[y*naux_2 +x]* eigval[x] * eigvec[z*naux_2+x];
        }


        libint2::Engine eri3_engine(libint2::Operator::coulomb,
                                    std::max(dfbs.max_nprim(),obs.max_nprim()),
                                    std::max(dfbs.max_l(),obs.max_l()));


        double* B = new double[naux_2*nroao*nroao];
        memset(B,0,naux_2*nroao*nroao*sizeof(double));

        eri3_engine.set_braket(libint2::BraKet::xs_xx);
        auto shell2bf_ob = obs.shell2bf();
        const auto& buf_vec_eri3 = eri3_engine.results();


        const auto& unitshell = libint2::Shell::unit();
        for(auto s1=0; s1!=dfbs.size(); ++s1) {
            for (auto s2=0; s2!= obs.size(); ++s2) {
                for (auto s3=0; s3!=obs.size(); ++s3) {
                    eri3_engine.compute(dfbs[s1],obs[s2],obs[s3]);
 //                   eri3_engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 0>(
//                                dfbs[s1], unitshell, obs[s2], obs[s3]);
                    auto ints_shellset_eri3 = buf_vec_eri3[0];

                    if (ints_shellset_eri3!=NULL) {
                        auto bf1 = shell2bf[s1]; // first basis function in first shell
                        auto n1 = dfbs[s1].size(); // number of basis functions in first shell

                        auto bf2 = shell2bf_ob[s2]; // first basis function in second shell
                        auto n2 =  obs[s2].size();// number of basis functions in second shell

                        auto bf3 = shell2bf_ob[s3]; // first basis function in second shell
                        auto n3 =  obs[s3].size();// number of basis functions in second shell

                        for(auto f1=0; f1!=n1; ++f1)
                            for(auto f2=0; f2!=n2; ++f2)
                                for(auto f3=0; f3!=n3; ++f3)
                                    B[(bf1+f1)*nroao*nroao + (bf2+f2)*nroao + (bf3+f3)] = ints_shellset_eri3[f1*n2*n3+f2*n3+f3];
                    }
                }
            }
        }
        BPQ  = new double[naux_2*nroao*nroao];
        double* BPQ2= new double[naux_2*nroao*nroao];
        memset(BPQ,0,naux_2*nroao*nroao*sizeof(double));


        //BPQ  = matmult(B,J^-1/2)
    //  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,naux_2,nroao*nroao,
     //               naux_2,1.0,J,naux_2,B,nroao*nroao,0.0,BPQ,nroao*nroao);

        for(int P=0; P<naux_2; P++) {
            for(int Q=0; Q<naux_2; Q++) {
                for(int pq=0; pq<nroao*nroao; pq++) {
                    BPQ[P*nroao*nroao + pq] += J[P*naux_2 + Q] * B[Q*nroao*nroao + pq];
                }
            }
        }


       // for (int i=0;i<naux_2*nroao*nroao;i++)
       //     diff+=BPQ2[i] - BPQ[i];

      //  std::cout<<"DIFF: "<<diff<<"\n";


        std::ifstream datf;
        datf.open(prefix+".rimp2");
        datf.read((char*) BPQ2, naux_2*nroao*nroao*sizeof(double));
        datf.close();

        double diff=0.0;
        for (int i=0;i<naux_2*nroao*nroao;i++)
            diff+=std::fabs(BPQ2[i] - BPQ[i]);

        std::cout<<"DIFF: "<<diff<<"\n";

        diff=0.0;
        for (int i=sortcount[0]; i<sortcount[1]; i++) {
            int a = intnums[4*i+0];
            int b = intnums[4*i+1];
            int c = intnums[4*i+2];
            int d = intnums[4*i+3];
            double integral = 0.0;
            for(int Q=0; Q<naux_2; Q++) {
                integral  += BPQ[Q*nroao*nroao+ a*nroao +b]*BPQ[Q*nroao*nroao+ c*nroao +d];
            }
            diff += std::fabs(integral - intval[i]);
            //std::cout<<a<<" "<<b<<" "<<" "<<c<<" "<<d<<" "<<integral<<" "<<intval[i]<<" "<< integral -intval[i] <<"\n";
        }
        std::cout<<"accumulated diffs:"<<diff;



        delete[] eigval;
        delete[] eigvec_c;
        delete[] eigvec;
        delete[] J;
        delete[] B;
    }
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

    double* FMo  = new double[nroao*nroao];
    build_FMo(nroao,Fmat,FMo,MOs);


    //RI-MP2 //trafo
    double ritrafo_start = 0.0;
    double ritrafo_end   = 0.0;
    double rimp2_start   = 0.0;
    double rimp2_end     = 0.0;
    {
        ritrafo_start = omp_get_wtime();
        double* Bia = new double[naux_2*nroe/2*(nroao-nroe/2)];
        read_transform_ri(prefix,nroe,nroao,naux_2,MOs,Bia);
        //transform_ri(nroe,nroao,naux_2,MOs,BPQ,Bia);
        ritrafo_end = omp_get_wtime();


        rimp2_start = omp_get_wtime();
        run_canonical_mp2_ri(nroe,nroao,naux_2,Bia,FMo);
        rimp2_end = omp_get_wtime();

        delete[] Bia;
    } //end RI-MP2

    //canonical mp2
    double mp2_start = omp_get_wtime();
    run_canonical_mp2(nroe,nroao,prec_ints,FMo);
    double mp2_end   = omp_get_wtime();


    std::cout <<"\n\nTIMINGS:\n";
    std::cout << std::setw( PWIDTH_L ) << "TEI [s]:" <<std::setw( PWIDTH_R )<<ints_end - ints_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "SCF [s]:" <<std::setw( PWIDTH_R )<<scf_end-scf_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "4-index [s]:" <<std::setw( PWIDTH_R )<<trafo_end-trafo_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "ri-4ind [s]:" <<std::setw( PWIDTH_R )<<ritrafo_end-ritrafo_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "RI-MP2 [s]:" <<std::setw( PWIDTH_R )<<rimp2_end-rimp2_start<< "s \n";
    std::cout << std::setw( PWIDTH_L ) << "MP2 [s]:" <<std::setw( PWIDTH_R )<<mp2_end-mp2_start<< "s \n";

    double total_end = omp_get_wtime();
    std::cout << "---------------------------" << '\n';
    std::cout<< std::setw( PWIDTH_L ) << "Total [s]:"<<std::setw( PWIDTH_R )<<total_end-total_start<<"s \n";
    std::cout << "\n\n--END--" << '\n';

    for(auto a: timer) {
        std::cout<<a.first<<" "<< a.second.second - a.second.first;
    }

    return 0;
}
