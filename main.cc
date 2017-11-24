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


	int nroe;                   //Nr of electrons  (if negative, read in center of mass)
	int llim;                   //first MO used for correlation
	int ulim;                //last  MO used for correlation
	int nroao;
	int nroa;
	long long int nrofint;    //Nr of two electron Integrals

	double* intval;
	unsigned short* intnums;
	long long int sortcount[4];

	double* dumd;
	double*        Dx;          //Dipole X                          nroao*nroao
	double*        Dy;            //Dipole Y                          nroao*nroao
	double*        Dz;            //Dipole Z                          nroao*nroao
	//atoms
	double*        coord;       //atomic coordinats               3*nroa
	double*        charges;     //atomic charges                    nroa
	double*        mass;        //atomic masses                     nroa

	//one electron mat&vecs
	double*        Smat;        //Overlap matrix S                  nroao*nroao
	double*        Hmat;        //one electron Hamiltionian         nroao*nroao
	double*        Tmat;        //Kinetic energy operator           nroao*nroao
	double*        Som12;       //S^-1/2                            nroao*nroao
	double*        Pmat_old;    //Old density matrix                nroao*nroao
	double*        Fmat;        //Fock matrix                       nroao*nroao
	double*        MOens;       //MO Energies
	double*        MOens_old;   //Old MO energies;                  nroao
	double* e1pmat;

	//Temporary memory spaces
	double*  tmpvecs;              //                               nroao*nroao
	double*  tmpvals;              //                               nroao
	double*  tmpmat1;              //                               nroao*nroao
	double*  tmpmat2;              //

	double mu_core[3];
	double*        Pmat;        //density matrix                    nroao*nroao
	double*        MOs;         //MO coeffs
	double ion_rep;

	double* prec_ints;


	std::string sysfile;
	std::string wavefile;

	double total_start = omp_get_wtime();

	print_header();
	if(argc != 4) {
		std::cerr << "Need PREFIX XYZ-file  Basis\n";
		exit(1);
	}

	sysfile  = std::string(argv[1])+".sys";
	wavefile = std::string(argv[1])+".ahfw";
	get_sys_size(sysfile,&nroe, &nroao, &nroa,  &nrofint);

	std::cout<< "\nSystem sizes read from       : " << sysfile << "\n";
	std::cout<< "Nr of basis functions        : " << nroao   << "\n";
	std::cout<< "Nr of electrons              : " << nroe   << "\n";
	std::cout<< "Nr of atoms                  : " << nroa << "\n";
	std::cout<< "Nr of non zero 2el integrals : " << nrofint << "\n";

	//MEMORY ALLOCATION for one electron atoms, mat&vecs
	int atom_ao_mem = 5*nroa+14*nroao*nroao+3*nroao;
	dumd  = (double *) malloc(atom_ao_mem*sizeof(double));
	int inc = 0;

	coord = &(dumd[inc]); inc += 3*nroa; charges = &(dumd[inc]); inc += nroa; mass  = &(dumd[inc]); inc += nroa;
	Smat = &(dumd[inc]); inc += nroao*nroao; Hmat = &(dumd[inc]); inc += nroao*nroao; Tmat = &(dumd[inc]); inc += nroao*nroao;
	Som12 = &(dumd[inc]); inc += nroao*nroao; Pmat = &(dumd[inc]); inc += nroao*nroao; Pmat_old = &(dumd[inc]); inc += nroao*nroao;
	Fmat = &(dumd[inc]); inc += nroao*nroao; MOs = &(dumd[inc]); inc += nroao*nroao; Dx = &(dumd[inc]); inc += nroao*nroao;
	Dy = &(dumd[inc]); inc += nroao*nroao; Dz = &(dumd[inc]); inc += nroao*nroao;

	MOens = &(dumd[inc]); inc += nroao; MOens_old = &(dumd[inc]); inc += nroao;
	tmpmat1 = &(dumd[inc]); inc += nroao*nroao; tmpmat2 = &(dumd[inc]); inc += nroao*nroao; tmpvecs = &(dumd[inc]); inc += nroao*nroao;

	tmpvals = &(dumd[inc]); inc += nroao;

	intval        = new double[nrofint];                   //two electron integrals
	intnums       = new unsigned short[nrofint*4];         //two electron indices
	//num of two electron Integrals in each perm. type
	read_sys(sysfile, coord, charges, mass, Hmat, Tmat, Smat, Dx, Dy,  Dz, sortcount, intval, intnums);

	std::cout<< "Atomic charges               : ";
	for (int i=0; i<nroa; i++)
		std::cout<<charges[i]<<" ";

	ion_rep =  calc_ion_rep( nroa, coord, charges);
	std::cout << "\nCalculating S^-1/2\n";
	std::cout << "Minimal eigenvalue of S is " << calc_S12( nroao, Smat, Som12, tmpmat1, tmpvecs, tmpvals) <<"\n";

	read_wav_HF(wavefile, nroao, MOens, MOs);
	std::cout<< "HF-wave function  read from " << wavefile << "\n";


	std::cout << "Starting libint2 - part" << '\n';
	libint2::initialize();

	std::string geometry_path(argv[2]);
	std::string basis(argv[3]);

	std::cout << "geometry file:" <<geometry_path<<  '\n';
	std::ifstream geometry_file(geometry_path);
	std::vector<libint2::Atom> atoms= libint2::read_dotxyz(geometry_file);
	libint2::BasisSet obs(basis,atoms);
	obs.set_pure(true);
	std::cout << "nbf:" <<obs.nbf()<<  '\n';


	libint2::Engine s_engine(libint2::Operator::overlap,obs.max_nprim(),obs.max_l());
	libint2::Engine t_engine(libint2::Operator::kinetic,obs.max_nprim(),obs.max_l());
	libint2::Engine v_engine(libint2::Operator::nuclear,obs.max_nprim(),obs.max_l());
	v_engine.set_params(libint2::make_point_charges(atoms));

//  std::copy(begin(obs), end(obs),
//          std::ostream_iterator<libint2::Shell>(std::cout, "\n"));


	double* Tmat_libint = new double[obs.nbf() * obs.nbf()];
	double* Vmat_libint = new double[obs.nbf() * obs.nbf()];
	double* Hmat_libint = new double[obs.nbf() * obs.nbf()];
	double* Smat_libint = new double[obs.nbf() * obs.nbf()];

	auto shell2bf = obs.shell2bf();
	std::vector<int> one_shift;
	std::vector<int> two_shift;
  int one_size = 0;
  int two_size = 0;
	const auto& buf_vec_t = t_engine.results();
	const auto& buf_vec_v = v_engine.results();
	const auto& buf_vec_s = s_engine.results();

	for(auto s1=0; s1!=obs.size(); ++s1) {
		for(auto s2=0; s2!=obs.size(); ++s2) {
			t_engine.compute(obs[s1], obs[s2]);
			v_engine.compute(obs[s1], obs[s2]);
      s_engine.compute(obs[s1], obs[s2]);

			auto ints_shellset_t = buf_vec_t[0];
			auto ints_shellset_v = buf_vec_v[0];
			auto ints_shellset_s = buf_vec_s[0];

			auto bf1 = shell2bf[s1]; // first basis function in first shell
			auto n1 = obs[s1].size(); // number of basis functions in first shell
			auto bf2 = shell2bf[s2]; // first basis function in second shell
			auto n2 = obs[s2].size(); // number of basis functions in second shell

			//libint order is -l,-l+1...0,1,...l
			if (obs[s1].contr[0].l == 0) {
				one_shift = {0};
        one_size  = 1;
			}
			if (obs[s2].contr[0].l == 0) {
				two_shift = {0};
        two_size  = 1;
			}

			if (obs[s1].contr[0].l == 1) {
				one_shift = {+2,-1, -1 };
        one_size  = 3;
			}
			if (obs[s2].contr[0].l == 1) {
				two_shift = {+2,-1, -1 };
        two_size  = 3;
			}
			if (obs[s1].contr[0].l == 2) {
				one_shift = {+4,+1,-2,-2,-1};
        one_size  = 5;
			}
			if (obs[s2].contr[0].l == 2) {
				two_shift = {+4,+1,-2,-2,-1};
        two_size  = 5;
			}

      if (obs[s1].contr[0].l == 3) {
				one_shift = {+6,+3,0,-3,-3,-2,-1};
        one_size  = 7;
			}
			if (obs[s2].contr[0].l == 3) {
				two_shift = {+6,+3,0,-3,-3,-2,-1};
        two_size  = 7;
			}

      if (obs[s1].contr[0].l == 4) {
				one_shift = {+8,+5,+2,-1,-4,-4,-3,-2,-1};
        one_size  = 9;
			}
			if (obs[s2].contr[0].l == 4) {
				two_shift = {+8,+5,+2,-1,-4,-4,-3,-2,-1};
        two_size  = 9;
			}
			// integrals are packed into ints_shellset in row-major (C) form
			// this iterates over integrals in this order
			for(auto f1=0; f1!=n1; ++f1)
				for(auto f2=0; f2!=n2; ++f2) {
					Tmat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_t[f1*n2+f2];
					Vmat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_v[f1*n2+f2];
					Smat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_s[f1*n2+f2];
				}
		}
	}


  double* Tmat_trans = new double[nroao*nroao];
  double* Smat_trans = new double[nroao*nroao];
  double* Hmat_trans = new double[nroao*nroao];

  for (size_t i = 0; i < nroao*nroao; i++)
    Tmat_trans[i] = 0;

  for(auto s1=0; s1!=obs.size(); ++s1) {
    for(auto s2=0; s2!=obs.size(); ++s2) {

      auto bf1 = shell2bf[s1]; // first basis function in first shell
			auto n1 = obs[s1].size(); // number of basis functions in first shell
			auto bf2 = shell2bf[s2]; // first basis function in second shell
			auto n2 = obs[s2].size(); // number of basis functions in second shell

      if (obs[s1].contr[0].l == 0) {
				one_shift = {0};
        one_size  = 1;
			}
			if (obs[s2].contr[0].l == 0) {
				two_shift = {0};
        two_size  = 1;
			}

			if (obs[s1].contr[0].l == 1) {
				one_shift = {+2,-1, -1 };
        one_size  = 3;
			}
			if (obs[s2].contr[0].l == 1) {
				two_shift = {+2,-1, -1 };
        two_size  = 3;
			}
			if (obs[s1].contr[0].l == 2) {
				one_shift = {+4,+1,-2,-2,-1};
        one_size  = 5;
			}
			if (obs[s2].contr[0].l == 2) {
				two_shift = {+4,+1,-2,-2,-1};
        two_size  = 5;
			}

      if (obs[s1].contr[0].l == 3) {
				one_shift = {+6,+3,0,-3,-3,-2,-1};
        one_size  = 7;
			}
			if (obs[s2].contr[0].l == 3) {
				two_shift = {+6,+3,0,-3,-3,-2,-1};
        two_size  = 7;
			}

      if (obs[s1].contr[0].l == 4) {
				one_shift = {+8,+5,+2,-1,-4,-4,-3,-2,-1};
        one_size  = 9;
			}
			if (obs[s2].contr[0].l == 4) {
				two_shift = {+8,+5,+2,-1,-4,-4,-3,-2,-1};
        two_size  = 9;
			}
      for(auto f1=0; f1!=n1; ++f1)
				for(auto f2=0; f2!=n2; ++f2) {
					Tmat_trans[(bf1+f1)*obs.nbf() + (bf2+f2)] = Tmat[(bf1+f1+one_shift[f1])*obs.nbf() + (bf2+f2+two_shift[f2])];
          Smat_trans[(bf1+f1)*obs.nbf() + (bf2+f2)] = Smat[(bf1+f1+one_shift[f1])*obs.nbf() + (bf2+f2+two_shift[f2])];
				}
		}
  }

    double tmatdiff = 0;
    double smatdiff = 0;
    double tmatmax = 0;
    double smatmax = 0;
    for (size_t i = 0; i < nroao*nroao; i++) {
      tmatdiff += fabs(Tmat_trans[i] - Tmat_libint[i]);
      smatdiff += fabs(Smat_trans[i] - Smat_libint[i]);
      if (fabs(Tmat_trans[i] - Tmat_libint[i])>tmatmax)
        tmatmax = fabs(Tmat_trans[i] - Tmat_libint[i]);
      if (fabs(Smat_trans[i] - Smat_libint[i]) > smatmax)
        smatmax = fabs(Smat_trans[i] - Smat_libint[i]);
    }


	std::cout << "tmatdiff:" <<tmatdiff << '\n';
	std::cout << "smatdiff:" <<smatdiff << '\n';
  std::cout << "tmatmax :" <<tmatmax << '\n';
	std::cout << "smatmax :" <<smatmax << '\n';
	//for (size_t i = 0; i < nroao; i++) {
	//  for (size_t j = 0; j < nroao; j++) {
	//  std::cout<<"I:"<<i<<" J:"<<j<<'\t'<<std::setw( PWIDTH_L )<<Tmat[i*nroao+j] <<std::setw( PWIDTH_L )<< Tmat_libint[i*nroao + j] << std::setw( PWIDTH_L )<<Vmat_libint[i*nroao + j] <<"\n";  /* code */
	//  }

	//}
	libint2::finalize();

//some core guess for testing
	// for (int i =0;i<nroao*nroao;i++)
	//	     Fmat[i] = Hmat[i];
//   double *tmpmat = new double[nroao*nroao];
	//diag_Fmat(nroao, Fmat,MOs,MOens,Som12, tmpmat);
	// delete[] tmpmat;
//  --END-COREGUESS
/*
   double testMOs = 0.0;
   for (int i=0;i<nroao;i++){
      testMOs += MOs[i*nroao+i];
   }
   if (testMOs == 0.0)
   {
      std::cout<<"\nForming core guess!\n";
      form_core_guess(nroao,Fmat,Hmat,Som12,MOs,MOens);
   }



   run_scf(nroao,nroe,MOs,Pmat,Hmat,Fmat,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);


   long long int prec_mem = (long long int) nroe/2*(long long int) nroao* (long long int) nroao* (long long int) nroao;
    std::cout << "Need " << prec_mem*sizeof(double) << " bytes (" << prec_mem*sizeof(double)/1024/1024 << " MB) for precalculation\n";
    std::cout << prec_mem << " doubles on "<< worldsize<< " cpus \n";
    if (worldsize > 1){
      std::cout << "Allocating " << prec_mem/worldsize << " on "<< (worldsize-1)<< " processors and "<< (prec_mem/worldsize + prec_mem%worldsize)<< " on one cpu";
    }

    if (rank==worldsize-1){
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
    for(i = 0; i < nroe/2; i++){
      for(j = 0; j < nroao; j++){
        for(k = 0; k < nroao; k++){
          for(l = 0; l <  nroao;l++){
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
    for (size_t i = 0; i < nroao*nroao; i++){
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

   double mp2_start = omp_get_wtime();
   run_canonical_mp2(nroe,nroao,prec_ints,FMo);
   double mp2_end   = omp_get_wtime();

   std::cout <<"\n\nTIMINGS:\n";
   std::cout << std::setw( PWIDTH_L ) << "4-index [s]:" <<std::setw( PWIDTH_R )<<trafo_end-trafo_start<< "s \n";
   std::cout << std::setw( PWIDTH_L ) << "MP2 [s]:" <<std::setw( PWIDTH_R )<<mp2_end-mp2_start<< "s \n";


   double total_end = omp_get_wtime();

   std::cout<< std::setw( PWIDTH_L ) << "Total [s]:"<<std::setw( PWIDTH_R )<<total_end-total_start<<"s \n";
   std::cout << "\n\n--END--" << '\n';
 */
	return 0;
}
