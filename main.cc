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

int main(int argc, char const *argv[]) {
  int worldsize = 1;
  int rank = 0;


  int           nroe;           //Nr of electrons  (if negative, read in center of mass)
  int           llim;           //first MO used for correlation
  int           ulim;        //last  MO used for correlation
  int           nroao;
  int           nroa;
  long long int nrofint;      //Nr of two electron Integrals

  double* intval;
	unsigned short* intnums;
  long long int  sortcount[4];

  double* dumd;
  double*        Dx;            //Dipole X                          nroao*nroao
	double*        Dy;            //Dipole Y                          nroao*nroao
	double*        Dz;            //Dipole Z                          nroao*nroao
  //atoms
  double*        coord;         //atomic coordinats               3*nroa
  double*        charges;       //atomic charges                    nroa
  double*        mass;          //atomic masses                     nroa

  //one electron mat&vecs
  double*        Smat;          //Overlap matrix S                  nroao*nroao
  double*        Hmat;          //one electron Hamiltionian         nroao*nroao
  double*        Tmat;          //Kinetic energy operator           nroao*nroao
  double*        Som12;         //S^-1/2                            nroao*nroao
  double*        Pmat_old;      //Old density matrix                nroao*nroao
  double*        Fmat;          //Fock matrix                       nroao*nroao
  double*        MOens;         //MO Energies
  double*        MOens_old;     //Old MO energies;                  nroao
  double* e1pmat ;

  //Temporary memory spaces
  double*  tmpvecs;                //                               nroao*nroao
  double*  tmpvals;                //                               nroao
  double*  tmpmat1;                //                               nroao*nroao
  double*  tmpmat2;                //

  double mu_core[3];
  double*        Pmat;          //density matrix                    nroao*nroao
  double*        MOs;           //MO coeffs
  double ion_rep;

  double* prec_ints;


  std::string sysfile;
  std::string wavefile;

  print_header();
  if(argc != 2){
    std::cerr << "Need prefix\n";
    exit(1);
  }
  sysfile  = std::string(argv[1])+".sys" ;
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
  Som12 = &(dumd[inc]); inc += nroao*nroao; Pmat = &(dumd[inc]); inc += nroao*nroao;Pmat_old = &(dumd[inc]); inc += nroao*nroao;
  Fmat = &(dumd[inc]); inc += nroao*nroao; MOs = &(dumd[inc]); inc += nroao*nroao; Dx = &(dumd[inc]); inc += nroao*nroao;
  Dy = &(dumd[inc]); inc += nroao*nroao; Dz = &(dumd[inc]); inc += nroao*nroao;

  MOens = &(dumd[inc]); inc += nroao; MOens_old = &(dumd[inc]); inc += nroao;
  tmpmat1 = &(dumd[inc]); inc += nroao*nroao; tmpmat2 = &(dumd[inc]); inc += nroao*nroao; tmpvecs = &(dumd[inc]); inc += nroao*nroao;

  tmpvals = &(dumd[inc]); inc += nroao;

  intval        = new double[nrofint];                     //two electron integrals
  intnums       = new unsigned short[nrofint*4];           //two electron indices
                                             //num of two electron Integrals in each perm. type
  read_sys(sysfile, coord, charges, mass, Hmat, Tmat, Smat, Dx, Dy,  Dz, sortcount, intval, intnums);

  std::cout<< "Atomic charges               : ";
  for (int i=0;i<nroa;i++)
	 std::cout<<charges[i]<<" ";

  ion_rep =  calc_ion_rep( nroa, coord, charges);
  std::cout << "\nCalculating S^-1/2\n";
  std::cout << "Minimal eigenvalue of S is " << calc_S12( nroao, Smat, Som12, tmpmat1, tmpvecs, tmpvals) <<"\n";

  read_wav_HF(wavefile, nroao, MOens, MOs);
  std::cout<< "HF-wave function  read from " << wavefile << "\n";

//some core guess for testing
 // for (int i =0;i<nroao*nroao;i++)
 //	     Fmat[i] = Hmat[i];
//   double *tmpmat = new double[nroao*nroao];
   //diag_Fmat(nroao, Fmat,MOs,MOens,Som12, tmpmat);
  // delete[] tmpmat;
//  --END-COREGUESS
  //form_core_guess(nroao,Fmat,Hmat,Som12,MOs,MOens);
  run_scf(nroao,nroe,MOs,Pmat,Hmat,Fmat,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);

  std::cout << "Precalculating <oo|oo> up to <ov|vv>\n";
  std::cout << nroe << '\n';
  std::cout << nroe/2 << '\n';
  std::cout << nroao << '\n';
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
    std::cout<<'\n'<<std::setw( 10 )<< "oneE:" <<std::setw( 16 )<<oneE<<'\n';

    for (size_t i = 0; i < nroe/2; i++) {
      for (size_t j = 0; j < nroe/2; j++) {
          twoE+= (2 *prec_ints[i*istep + i*jstep + j*kstep + j]
                   - prec_ints[i*istep + j*jstep + i*kstep + j]);

      }
    }
    std::cout << std::setw( 10 ) << "twoE:" <<std::setw( 16 )<<twoE<< '\n';
    std::cout << std::setw( 10 ) << "HF  :" <<std::setw( 16 )<<oneE+twoE<<'\n';
    ///END Hartree-Fock

  run_canonical_mp2(nroe,nroao,prec_ints,FMo);

  std::cout << "\n\n--END--" << '\n';
  return 0;
}
