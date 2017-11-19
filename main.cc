#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "IO/ops_io.h"
#include "SCF/ops_rhf.h"
#include "SCF/scf.h"
#include "MOTRANS/ops_cis.h"
#include <iomanip>

void read_input(std::ifstream* inputfile, std::string* sysfile, int *nroe,int *llim, int *ulim, std::string* wavefile)
{
  std::string tmpline;
  std::stringstream ss_input;
  if (inputfile->is_open())
  {
    while(std::getline((*inputfile),tmpline))
    {
      //remove comments
      if (tmpline.find_first_of("#")!=std::string::npos)
        tmpline.erase(tmpline.find_first_of("#"),tmpline.length());
      if (tmpline.find_first_of(" ")!=std::string::npos)
        tmpline.erase(tmpline.find_first_of(" "),tmpline.length());
      if (tmpline.compare("\n")!=0)
        ss_input <<tmpline<<" ";
    }
  }
  inputfile->close();
  ss_input >> (*sysfile) >> (*nroe) >> (*llim) >> (*ulim) >> (*wavefile);
  std::cout<<"Sysfile    : "<<(*sysfile)<<std::endl;
  std::cout<<"#Electrons : "<<(*nroe)<<std::endl;
  std::cout<<"llim       : "<<(*llim)<<std::endl;
  std::cout<<"ulim       : "<<(*ulim)<<std::endl;
  std::cout<<"Wavefile   : "<<(*wavefile)<<std::endl;
  std::cout<<"=======================================\n\n";
}


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

  //auto mymol  = psi::molecule::

  std::string sysfile;
  std::string wavefile;

  std::cout<<"PSI4-MPI Program \n";
  if(argc != 3){
    std::cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  std::ifstream inputfile (argv[1]);

  read_input(&inputfile,&sysfile,&nroe,&llim,&ulim,&wavefile);
  get_sys_size(sysfile, &nroao, &nroa,  &nrofint);

  std::cout<< "System sizes read from       : " << sysfile << "\n";
  std::cout<< "Nr of basis functions        : " << nroao   << "\n";
  std::cout<< "Nr of atoms                  : " << nroa << "\n";
  std::cout<< "Nr of non zero 2el integrals : " << nrofint << "\n";
  std::cout<< "=======================================\n\n";

  //MEMORY ALLOCATION for one electron atoms, mat&vecs
  int atom_ao_mem = 5*nroa+14*nroao*nroao+3*nroao;
  std::cout << "Need " << atom_ao_mem*sizeof(double) << " bytes for atomic + one electron data\n";

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
  //MEMORY ALLOCATION for two electron values
  std::cout << "Need " << nrofint*(sizeof(double)+sizeof(unsigned short)*4) << " bytes for two electron data\n";

  intval        = new double[nrofint];                     //two electron integrals
  intnums       = new unsigned short[nrofint*4];           //two electron indices
                                             //num of two electron Integrals in each perm. type
  read_sys(sysfile, coord, charges, mass, Hmat, Tmat, Smat, Dx, Dy,  Dz, sortcount, intval, intnums);
  ion_rep =  calc_ion_rep( nroa, coord, charges);

  std::cout << "Calculating S^-1/2\n";
  std::cout << "Minimal eigenvalue of S is " << calc_S12( nroao, Smat, Som12, tmpmat1, tmpvecs, tmpvals) <<"\n";

  read_wav_HF(wavefile, nroao, MOens, MOs);
  std::cout<< "HF-wave function  read from " << wavefile << "\n";

  for (size_t x = 0; x < nroao*nroao; x++) {
    Fmat[x] = Hmat[x];
  }
  run_scf(nroao,nroe,MOs,Pmat,Hmat,Fmat,intnums,intval,sortcount,nrofint,Som12,100,ion_rep);





    //
    //  TRAFO
    //
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

/*
    std::cout << "\nCalculating semi-canonical amplitudes ... ";
    int nocc = nroe/2;
    int nvir = nroao - nocc;
    double* T_ijab   = new double[nocc*nocc*nvir*nvir];
    int Ti,Tj,Ta,Tb = 0;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    for (Ti = 0; Ti < nocc; Ti++) {
      for (Tj = 0; Tj < nocc; Tj++) {
        for (Ta = 0; Ta < nvir; Ta++) {
          for (Tb = 0; Tb < nvir; Tb++) {
            T_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] = - prec_ints[Ti * istep + (Ta+nocc)*jstep  + Tj*kstep + (Tb+nocc)] /
                                                                  (FMo[(nocc+Ta)*nroao + (nocc+Ta)] + FMo[(nocc+Tb)*nroao + (nocc+Tb)] - FMo[Ti*nroao + Ti] - FMo[Tj*nroao + Tj]);

          }
        }
      }
    }
    std::cout << "done" << '\n';
    double* G_ijab   = new double[nocc*nocc*nvir*nvir];
    double* R_ijab   = new double[nocc*nocc*nvir*nvir];

    for (size_t iter = 0; iter < 20; iter++) {
      for (Ti = 0; Ti < nocc; Ti++) {
        for (Tj = 0; Tj < nocc; Tj++) {
          for (Ta = 0; Ta < nvir; Ta++) {
            for (Tb = 0; Tb < nvir; Tb++) {
              G_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] = 0.0;
              for (int k = 0; k < nocc; k++) {
                if (k!=Ti){
                  G_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] +=  -FMo[Ti*nroao + k] * T_ijab[k*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb];
                }
              }
            }
          }
        }
      }


      double Rsum = 0.0;
      for (Ti = 0; Ti < nocc; Ti++) {
        for (Tj = 0; Tj < nocc; Tj++) {
          for (Ta = 0; Ta < nvir; Ta++) {
            for (Tb = 0; Tb < nvir; Tb++) {
              R_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] = prec_ints[Ti * istep + (Ta+nocc)*jstep  + Tj*kstep + (Tb+nocc)] +
                                                                (FMo[(nocc+Ta)*nroao + (nocc+Ta)] + FMo[(nocc+Tb)*nroao + (nocc+Tb)] - FMo[Ti*nroao + Ti] - FMo[Tj*nroao + Tj])*T_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb]+
                                                                G_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] +
                                                                G_ijab[Tj*Ti_step + Ti*Tj_step + Tb*Ta_step + Ta];
              Rsum += fabs(R_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb]);
          }
        }
      }
    }
    double E = 0.0;
    for (Ti = 0; Ti < nocc; Ti++) {
      for (Tj = 0; Tj < nocc; Tj++) {
        for (Ta = 0; Ta < nvir; Ta++) {
          for (Tb = 0; Tb < nvir; Tb++) {
            E +=  (prec_ints[Ti * istep + (Ta+nocc)*jstep  + Tj*kstep + (Tb+nocc)] + R_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb])* (2 * T_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] -T_ijab[Ti*Ti_step + Tj*Tj_step + Tb*Ta_step + Ta]);
          }
        }
      }
    }

    for (Ti = 0; Ti < nocc; Ti++) {
      for (Tj = 0; Tj < nocc; Tj++) {
        for (Ta = 0; Ta < nvir; Ta++) {
          for (Tb = 0; Tb < nvir; Tb++) {
              T_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb] += - R_ijab[Ti*Ti_step + Tj*Tj_step + Ta*Ta_step + Tb]/(FMo[(nocc+Ta)*nroao + (nocc+Ta)] + FMo[(nocc+Tb)*nroao + (nocc+Tb)] - FMo[Ti*nroao + Ti] - FMo[Tj*nroao + Tj]);
          }
        }
      }
    }
    std::cout <<std::fixed<<std::setw( 10 )<<std::setprecision(10)<<"E(sem-loc):" <<std::setw( 16 ) <<E<<'\t' << "Rsum:" << Rsum <<'\n';
  }
*/
    //Hartree-Fock
    for (size_t i = 0; i < nroe/2; i++) {
      oneE += 2 * HMo[i*nroao + i];
    }
    oneE += ion_rep;
    std::cout<<std::setw( 10 )<< "oneE:" <<std::setw( 16 )<<oneE<<'\n';

    for (size_t i = 0; i < nroe/2; i++) {
      for (size_t j = 0; j < nroe/2; j++) {
          twoE+= (2 *prec_ints[i*istep + i*jstep + j*kstep + j]
                   - prec_ints[i*istep + j*jstep + i*kstep + j]);

      }
    }
    std::cout << std::setw( 10 ) << "twoE:" <<std::setw( 16 )<<twoE<< '\n';
    std::cout << std::setw( 10 ) << "HF  :" <<std::setw( 16 )<<oneE+twoE<<'\n';
    ///END Hartree-Fock
/*
    //calculate canonical MP2:
    double EMP2 = 0.0;
    double EMP2_SS = 0.0;
    double EMP2_OS = 0.0;
    double* eps = new double[nroao];
    for (int i = 0; i < nroao; i++) {
      eps[i] = HMo[i*nroao + i];

      for (size_t j = 0; j < nroe/2; j++) {
        eps[i]+=(2 *prec_ints[j*istep + j*jstep + i*kstep + i]
                 - prec_ints[j*istep + i*jstep + j*kstep + i]);
      }
    }

    for (int i = 0; i < nroe/2; i++) {
      for (int j = 0; j < nroe/2; j++) {
        for (int a = nroe/2; a < nroao; a++) {
          for (int b = nroe/2;  b < nroao; b++) {
            EMP2_SS += - (prec_ints[i*istep + a*jstep + j*kstep + b]*(prec_ints[i*istep + a*jstep + j*kstep + b]))/ (eps[a]+eps[b]-eps[i]-eps[j]);
            EMP2_OS += - (prec_ints[i*istep + a*jstep + j*kstep + b] - prec_ints[i*istep + b*jstep + j*kstep + a])*prec_ints[i*istep + a*jstep + j*kstep + b]/(eps[a]+eps[b]-eps[i]-eps[j]);
          }
        }
      }
    }
    EMP2 = EMP2_OS + EMP2_SS;
    std::cout <<std::setw( 10 ) << "EMP2_SS:" <<std::setw( 16 )<<EMP2_SS<< '\n';
    std::cout <<std::setw( 10 ) << "EMP2_OS:" <<std::setw( 16 )<<EMP2_OS<< '\n';
    std::cout <<std::setw( 10 ) << "EMP2:"    <<std::setw( 16 )<<EMP2<< '\n';
*/








  std::cout << "\n\n--END--" << '\n';
  return 0;
}
