/*
 *@BEGIN LICENSE
 *
 * syshfwplug by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/factory.h"
#include "psi4/lib3index/cholesky.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/psimath.h"
#include "psi4/libqt/qt.h"
#include "psi4/lib3index/dftensor.h"
#include <fstream>
#include <iomanip>

void resort_integrals(unsigned short int*  intnums, double* intval,  long long int nrofint, long long int* sortcount );

void print_header(){
	std::cout << "++++++++++++++++++++++ " << '\n';
	std::cout << "+SYSHFWPLUG - plugin + " << '\n';
	std::cout << "++++++++++++++++++++++ " << '\n';
}






namespace psi {

namespace syshfwplug  {

void run_dryrun(SharedWavefunction ref_wfn, bool do_tei, bool do_rijk, bool do_rimp2){
	std::cout << "\nDRYRUN!:" << '\n';
	std::cout << "--------"<<"\n";

	auto molecule = ref_wfn->molecule();
	auto aoBasis  = ref_wfn->basisset();
	auto aux      = ref_wfn->get_basisset("JKFIT");
	auto aux2     = ref_wfn->get_basisset("RIFIT");

	int nbf   = aoBasis->nbf();
	int naux  = aux->nbf();
	int naux2 = aux2->nbf();


	std::cout << "NBFS" << '\n';
	std::cout << "AO    :" <<nbf<< '\n';
	std::cout << "RIJK  :" <<naux<< '\n';
	std::cout << "RIMP2 :" <<naux2<< '\n';


	if (do_tei) {

		std::shared_ptr<IntegralFactory> integral(new IntegralFactory(aoBasis, aoBasis, aoBasis, aoBasis));
		const double cutoff2el = 1.e-12;
		long long int count=0;

		std::shared_ptr<TwoBodyAOInt> eri(integral->eri());
		const double *buffer = eri->buffer();
		AOShellCombinationsIterator shellIter = integral->shells_iterator();
		for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
			// Compute quartet
			eri->compute_shell(shellIter);
			// From the quartet get all the integrals
			AOIntegralsIterator intIter = shellIter.integrals_iterator();
			for (intIter.first(); intIter.is_done() == false; intIter.next()) {
				int p = intIter.i();
				int q = intIter.j();
				int r = intIter.k();
				int s = intIter.l();
				if( fabs(buffer[intIter.index()]) > cutoff2el) {
					++count;
					if (count%(long long int)1.25E8 == 0) {
						std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
					}
				}
			}
		}
		std::cout<<"\nORBBAS\nSieve: "<<count*8<<" bytes\n"<<"Total: "<<(long long int) nbf*nbf*nbf*nbf*8<<" bytes";
		std::cout<<"\nSieve: "<<count/1024/1024/1024*8<<" GB \n"<<"Total: "<<(long long int) nbf*nbf*nbf*nbf*8/1024/1024/1024<<" GB \n\n";
	}

	if (do_rijk) {
		const double cutoff2el = 1.e-12;
		long long int count=0;
		auto zero = BasisSet::zero_ao_basis_set();
		std::shared_ptr<IntegralFactory> fact(new IntegralFactory(aux,zero,aoBasis,aoBasis));
		std::shared_ptr<TwoBodyAOInt> auxeri(fact->eri());
		const double *buffer = auxeri->buffer();
		count = 0;
		for (int P = 0; P < aux->nshell(); P++) {
			int np = aux->shell(P).nfunction();
			int pstart = aux->shell(P).function_index();
			for (int M = 0; M < aoBasis->nshell(); M++) {
				int nm = aoBasis->shell(M).nfunction();
				int mstart = aoBasis->shell(M).function_index();
				for (int N = 0; N < aoBasis->nshell(); N++) {
					int nn = aoBasis->shell(N).nfunction();
					int nstart = aoBasis->shell(N).function_index();

					auxeri->compute_shell(P,0,M,N);
					for (int p = 0, index = 0; p < np; p++) {
						for (int m = 0; m < nm; m++) {
							for (int n = 0; n < nn; n++, index++) {
								//Bp[p + pstart][(m + mstart) * nbf + (n + nstart)] = buffer[index];
								if (fabs(buffer[index])  > cutoff2el) {
									count++;
									if (count%(long long int)1.25E8 == 0) {
										std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
									}
								}
							}
						}
					}
				}
			}
		}
		std::cout<<"JKFIT\nSieve: "<<count*8<<" bytes\n"<<"Total: "<<(long long int)nbf*nbf*naux*8<<" bytes";
		std::cout<<"\nSieve: "<<count*8/1024/1024/1024<<" GB\n"<<"Total: "<<(long long int)nbf*nbf*naux*8/1024/1024/1024<<" GB\n\n";

	}

	if (do_rimp2) {
		const double cutoff2el = 1.e-12;
		long long int count=0;
		auto zero = BasisSet::zero_ao_basis_set();
		std::shared_ptr<IntegralFactory> fact(new IntegralFactory(aux2,zero,aoBasis,aoBasis));
		std::shared_ptr<TwoBodyAOInt> auxeri(fact->eri());
		const double *buffer = auxeri->buffer();
		count = 0;
		for (int P = 0; P < aux2->nshell(); P++) {
			int np = aux2->shell(P).nfunction();
			int pstart = aux2->shell(P).function_index();
			for (int M = 0; M < aoBasis->nshell(); M++) {
				int nm = aoBasis->shell(M).nfunction();
				int mstart = aoBasis->shell(M).function_index();
				for (int N = 0; N < aoBasis->nshell(); N++) {
					int nn = aoBasis->shell(N).nfunction();
					int nstart = aoBasis->shell(N).function_index();

					auxeri->compute_shell(P,0,M,N);
					for (int p = 0, index = 0; p < np; p++) {
						for (int m = 0; m < nm; m++) {
							for (int n = 0; n < nn; n++, index++) {
								//Bp[p + pstart][(m + mstart) * nbf + (n + nstart)] = buffer[index];
								if (fabs(buffer[index])  > cutoff2el) {
									count++;
									if (count%(long long int)1.25E8 == 0) {
										std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
									}
								}
							}
						}
					}
				}
			}
		}
		std::cout<<"RIMP2\nSieve: "<<count*8<<" bytes\n"<<"Total: "<<(long long int)nbf*nbf*naux2*8<<" bytes";
		std::cout<<"\nSieve: "<<count*8/1024/1024/1024<<" GB\n"<<"Total: "<<(long long int)nbf*nbf*naux2*8/1024/1024/1024<<" GB\n\n";
	}
}

void run_oei(SharedWavefunction ref_wfn, std::string pref){
	std::cout << "\nONE-ELECTRON INTEGRALS" << '\n';
	std::cout << "----------------------" << '\n';


	auto molecule  = ref_wfn->molecule();
	auto aoBasis   = ref_wfn->basisset();
	int nbf = aoBasis->nbf();


	std::shared_ptr<IntegralFactory> integral(new IntegralFactory(aoBasis, aoBasis, aoBasis, aoBasis));
	std::shared_ptr<MatrixFactory> factory(new MatrixFactory);
	factory->init_with(1, &nbf, &nbf);

	std::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
	std::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
	std::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
	std::shared_ptr<OneBodyAOInt> ecpOBI(integral->ao_ecp());
	std::shared_ptr<OneBodyAOInt> dip(integral->ao_dipole());

	std::shared_ptr<Matrix> sMat(factory->create_matrix("Overlap"));
	std::shared_ptr<Matrix> tMat(factory->create_matrix("Kinetic"));
	std::shared_ptr<Matrix> vMat(factory->create_matrix("Potential"));
	std::shared_ptr<Matrix> hMat(factory->create_matrix("One Electron Ints"));
	std::shared_ptr<Matrix> ecpMat(factory->create_matrix("EPC Matrix "));

	sOBI->compute(sMat);

	if (ref_wfn->reference_energy() != 0.0)
	{
		std::cout << "Reference Found!" << '\n';
		hMat   = ref_wfn->H();
		tOBI->compute(tMat);
		vOBI->compute(vMat);
	}
	else
	{
		std::cout << "No Reference Found! Calculating the OEIs" << '\n';
		tOBI->compute(tMat);
		vOBI->compute(vMat);
		ecpOBI->compute(ecpMat);
		hMat->copy(tMat);
		hMat->add(vMat);
		hMat->add(ecpMat);
	}

	//prepare output
	int nroao = nbf;

	double* hmat = new double[nroao*nroao];
	double* kmat = new double[nroao*nroao];
	double* smat = new double[nroao*nroao];
	double* vmat = new double[nroao*nroao];
	for(int x = 0; x < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			hmat[x*nroao+y] = hMat->get(x,y);
			kmat[x*nroao+y] = tMat->get(x,y);
			smat[x*nroao+y] = sMat->get(x,y);
			vmat[x*nroao+y] = vMat->get(x,y);
		}
	}

	std::ofstream datf;
	datf.open(pref+".oei");
	datf.write((char * ) hmat, sizeof(double)*nroao*nroao);
	datf.write((char * ) kmat, sizeof(double)*nroao*nroao);
	datf.write((char * ) smat, sizeof(double)*nroao*nroao);
	datf.write((char * ) vmat, sizeof(double)*nroao*nroao);
	datf.close();
	delete[] hmat;
	delete[] kmat;
	delete[] smat;
	delete[] vmat;

	std::cout << "OEI successfully exported" << '\n';
}


long long int run_tei(SharedWavefunction ref_wfn, std::string pref) {
	std::cout << "\nTEI-ELECTRON INTEGRALS" << '\n';
	std::cout << "----------------------" << '\n';

	auto molecule = ref_wfn->molecule();
	auto aoBasis  = ref_wfn->basisset();
	const double cutoff2el = 1.e-12;
	long long int count=0;

	std::shared_ptr<IntegralFactory> integral(new IntegralFactory(aoBasis, aoBasis, aoBasis, aoBasis));


	std::shared_ptr<TwoBodyAOInt> eri(integral->eri());
	const double *buffer = eri->buffer();
	AOShellCombinationsIterator shellIter = integral->shells_iterator();
	for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
		// Compute quartet
		eri->compute_shell(shellIter);
		// From the quartet get all the integrals
		AOIntegralsIterator intIter = shellIter.integrals_iterator();
		for (intIter.first(); intIter.is_done() == false; intIter.next()) {
			int p = intIter.i();
			int q = intIter.j();
			int r = intIter.k();
			int s = intIter.l();
			if( fabs(buffer[intIter.index()]) > cutoff2el) {
				++count;
				if (count%(long long int)1.25E8 == 0) {
					std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
				}
			}
		}
	}
	long long int nrofint = count;
	long long int sortcount[4] = {0, 0, 0, 0};
	double* intval = new double[nrofint];
	unsigned short* intnums = new unsigned short[nrofint*4+4];
	count = 0;
	for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
		// Compute quartet
		eri->compute_shell(shellIter);
		// From the quartet get all the integrals
		AOIntegralsIterator intIter = shellIter.integrals_iterator();
		for (intIter.first(); intIter.is_done() == false; intIter.next()) {
			int p = intIter.i();
			int q = intIter.j();
			int r = intIter.k();
			int s = intIter.l();
			if( fabs(buffer[intIter.index()]) > cutoff2el) {
				//	  outfile->Printf("\t(%2d %2d | %2d %2d) = %20.15f\n",
				//			  p, q, r, s, buffer[intIter.index()]);
				intval[count] = buffer[intIter.index()];
				intnums[4*count+0] = p;
				intnums[4*count+1] = q;
				intnums[4*count+2] = r;
				intnums[4*count+3] = s;
				++count;
			}
		}
	}
	resort_integrals(intnums, intval,  nrofint, sortcount);
	std::ofstream datf;
	datf.open(pref+".tei");
	datf.write((char *) sortcount, sizeof(long long int)*4);
	datf.write((char *) intval,    sizeof(double)*nrofint);
	datf.write((char *) intnums,   sizeof(unsigned short)*nrofint*4);
	datf.close();

	std::cout << "TEI successfully written" << '\n';
	return nrofint;
}


long long int run_rijk(SharedWavefunction ref_wfn, std::string pref){
	std::cout << "\nRIJK INTEGRALS:" << '\n';
	std::cout <<   "---------------"<<"\n";

	auto molecule = ref_wfn->molecule();
	auto aoBasis  = ref_wfn->basisset();
	auto aux      = ref_wfn->get_basisset("JKFIT");

	const double cutoff2el = 1.e-12;
	long long int count=0;
	auto zero = BasisSet::zero_ao_basis_set();
	std::shared_ptr<IntegralFactory> fact(new IntegralFactory(aux,zero,aoBasis,aoBasis));
	std::shared_ptr<TwoBodyAOInt> auxeri(fact->eri());
	const double *buffer = auxeri->buffer();
	count = 0;
	for (int P = 0; P < aux->nshell(); P++) {
		int np = aux->shell(P).nfunction();
		int pstart = aux->shell(P).function_index();
		for (int M = 0; M < aoBasis->nshell(); M++) {
			int nm = aoBasis->shell(M).nfunction();
			int mstart = aoBasis->shell(M).function_index();
			for (int N = 0; N < aoBasis->nshell(); N++) {
				int nn = aoBasis->shell(N).nfunction();
				int nstart = aoBasis->shell(N).function_index();

				auxeri->compute_shell(P,0,M,N);
				for (int p = 0, index = 0; p < np; p++) {
					for (int m = 0; m < nm; m++) {
						for (int n = 0; n < nn; n++, index++) {
							//Bp[p + pstart][(m + mstart) * nbf + (n + nstart)] = buffer[index];
							if (fabs(buffer[index])  > cutoff2el) {
								count++;
								if (count%(long long int)1.25E8 == 0) {
									std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
								}
							}
						}
					}
				}
			}
		}
	}

	std::cout<<"Export not yet implemented\n";
	return count;
}


long long int run_rimp2(SharedWavefunction ref_wfn, std::string pref){
	std::cout << "\nRIMP2 INTEGRALS:" << '\n';
	std::cout <<   "---------------"<<"\n";


	auto aoBasis  = ref_wfn->basisset();
	auto aux2     = ref_wfn->get_basisset("RIFIT");

	int nbf   = aoBasis->nbf();
	int naux2 = aux2->nbf();

	const double cutoff2el = 1.e-12;
	long long int count=0;
	auto zero = BasisSet::zero_ao_basis_set();
	std::shared_ptr<IntegralFactory> fact(new IntegralFactory(aux2,zero,aoBasis,aoBasis));
	std::shared_ptr<TwoBodyAOInt> auxeri(fact->eri());
	const double *buffer = auxeri->buffer();
	count = 0;
	for (int P = 0; P < aux2->nshell(); P++) {
		int np = aux2->shell(P).nfunction();
		int pstart = aux2->shell(P).function_index();
		for (int M = 0; M < aoBasis->nshell(); M++) {
			int nm = aoBasis->shell(M).nfunction();
			int mstart = aoBasis->shell(M).function_index();
			for (int N = 0; N < aoBasis->nshell(); N++) {
				int nn = aoBasis->shell(N).nfunction();
				int nstart = aoBasis->shell(N).function_index();

				auxeri->compute_shell(P,0,M,N);
				for (int p = 0, index = 0; p < np; p++) {
					for (int m = 0; m < nm; m++) {
						for (int n = 0; n < nn; n++, index++) {
							//Bp[p + pstart][(m + mstart) * nbf + (n + nstart)] = buffer[index];
							if (fabs(buffer[index])  > cutoff2el) {
								count++;
								if (count%(long long int)1.25E8 == 0) {
									std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
								}
							}
						}
					}
				}
			}
		}
	}

	std::cout << "RIMP2 export not yet implemented" << '\n';
	return count;



}


void run_export_mos(SharedWavefunction ref_wfn,std::string pref){
	int nroao = ref_wfn->basisset()->nbf();
	auto Calpha = std::make_shared<Matrix>(nroao,nroao);
	auto Cbeta = std::make_shared<Matrix>(nroao,nroao);


	double* MOensA  = new double[nroao];   //nroao
	double* MOensB  = new double[nroao];   //nroao
	double* MOsA    = new double[nroao*nroao];   //nroao*nroao
	double* MOsB    = new double[nroao*nroao]; //nroao*nroao

	if (ref_wfn->reference_energy() != 0.0)
	{
		std::cout << "Reference Found!" << '\n';
		Calpha = ref_wfn->Ca();
		Cbeta  = ref_wfn->Cb();
	}
	//maybe later
	for(int x = 0; x < nroao; x++)
		MOensA[x] = MOensB[x] = 0.;

	for(int x = 0; x < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			MOsA[x*nroao+y] = Calpha->get(y,x);
			MOsB[x*nroao+y] = Cbeta->get(y,x);
		}
	}

	std::ofstream datf;
	datf.open(pref+".ahfw");
	datf.write((char *)  &nroao, sizeof(int));
	datf.write((char * ) MOensA, sizeof(double)*nroao);
	datf.write((char * ) MOsA, sizeof(double)*nroao*nroao);
	datf.close();
	std::cout << "alpha Orbitals successfully exported" << '\n';

	datf.open(pref+".bhfw");
	datf.write((char *)  &nroao, sizeof(int));
	datf.write((char * ) MOensB, sizeof(double)*nroao);
	datf.write((char * ) MOsB, sizeof(double)*nroao*nroao);
	datf.close();
	std::cout << "beta Orbitals successfully exported" << '\n';
	delete[] MOsA;
	delete[] MOsB;
	delete[] MOensA;
	delete[] MOensB;
}

void run_export_sys(SharedWavefunction ref_wfn,std::string pref,
										long long int nrofint,
										long long int nrofaux,
										long long int nrofaux2){
	int nroe    = ref_wfn->nalpha() + ref_wfn->nbeta();
	int nroao   = ref_wfn->basisset()->nbf();
	int nroa  = ref_wfn->molecule()->natom();
	int naux_1  = ref_wfn->get_basisset("JKFIT")->nbf();
	int naux_2  = ref_wfn->get_basisset("RIFIT")->nbf();
	auto molecule = ref_wfn->molecule();
	std::shared_ptr<Matrix> coord (new Matrix(molecule->geometry()));
	std::string basisNameOB = ref_wfn->basisset()->name();
	std::string basisNameJK = ref_wfn->get_basisset("JKFIT")->name();
	std::string basisNameRI = ref_wfn->get_basisset("RIFIT")->name();

	basisNameOB.resize(32);
	basisNameJK.resize(32);
	basisNameRI.resize(32);



	double* Coord   = new double[3*nroa];
	int* mass       = new int[nroa];
	double* charges = new double[nroa];


	for (int j=0; j<molecule->natom(); j++) {
		charges[j] = molecule->Z(j);
		mass[j]    = molecule->mass(j);
		for (int i=0; i<3; i++) {
			Coord[3*j+i] = coord->get(j,i);
		}
	}

	std::cout << "\nSYSTEMDATA" << '\n';
	std::cout << "----------" << "\n\n";
	for (size_t i = 0; i < nroa; i++) {
    std::cout<<std::setw(2)<<charges[i]<<std::setw(10)<< Coord[i*3+0]<<std::setw(10) << Coord[i*3+1]<<std::setw(10)<< Coord[i*3+2]<<'\n';
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

	std::ofstream datf;

	datf.open(pref+".sys");
	datf.write((char *) &nroe, sizeof(int));
	datf.write((char *) &nroa, sizeof(int));
	datf.write((char *) &nroao, sizeof(int));
	datf.write((char *) &naux_1, sizeof(int));
	datf.write((char *) &naux_2, sizeof(int));
	datf.write((char *) &nrofint,   sizeof(long long int));
	datf.write((char *) &nrofaux,   sizeof(long long int));
	datf.write((char *) &nrofaux2,  sizeof(long long int));
	datf.write((char *) Coord, sizeof(double)*3*nroa);
	datf.write((char *) charges, sizeof(double)*nroa);
	datf.write((char *) mass,    sizeof(double)*nroa);
	datf.write((char *) basisNameOB.c_str(),basisNameOB.size());
	datf.write((char *) basisNameJK.c_str(),basisNameJK.size());
	datf.write((char *) basisNameRI.c_str(),basisNameRI.size());
	datf.close();
	std::cout << "Systeminfo sucessfully written" << '\n';
}

extern "C"
int read_options(std::string name, Options &options)
{
	if (name == "SYSHFWPLUG"|| options.read_globals()) {
		/*- The amount of information printed
		    to the output file -*/
		options.add_int("PRINT", 1);
		/*- Whether to compute two-electron integrals -*/
		options.add_bool("DO_TEI", true);
		options.add_bool("DO_OEI", true);
		options.add_bool("DRYRUN", false);
		options.add_bool("DO_RIJK", false);
		options.add_bool("DO_RIMP2", false);
		options.add_str("PREF","vergessen-wa");

	}

	return true;
}

extern "C"
SharedWavefunction syshfwplug(SharedWavefunction ref_wfn, Options &options)
{
	// Grab options from the options object
	int print         = options.get_int("PRINT");
	std::string pref  = options.get_str("PREF");
	bool dryrun       = options.get_bool("DRYRUN");
	bool do_rijk      = options.get_bool("DO_RIJK");
	bool do_rimp2     = options.get_bool("DO_RIMP2");
	bool do_tei       = options.get_bool("DO_TEI");
	bool do_oei       = options.get_bool("DO_OEI");

	print_header();
	std::cout << "Options:" << '\n';
	std::cout << "--------" << '\n';
	std::cout << "DO_OEI      :" <<do_oei<< '\n';
	std::cout << "DO_TEI      :" <<do_tei<< '\n';
	std::cout << "RIJK        :" <<do_rijk<< '\n';
	std::cout << "RIMP2       :" <<do_rimp2<< '\n';
	std::cout << "DRYRUN      :" <<dryrun<< '\n';
	std::cout << "PREF        :" <<pref<< '\n';


	if (dryrun)
	{
		run_dryrun(ref_wfn,do_tei,do_rijk,do_rimp2);
		return ref_wfn;
	}

	std::cout << "\nPRODUCTIVE RUN:" << '\n';
	std::cout << "---------------"<<"\n";

	long long int nrofint  = 0;              //nr of nonzero to electron integrals
	long long int nrofaux  = 0;              //nr of nonzero rijk integrals
	long long int nrofaux2 = 0;              //nr of nonzero rimp2 integrals




	if (do_oei) {
		run_oei(ref_wfn,pref);
	}

	if (do_tei)
	{
		nrofint = run_tei(ref_wfn,pref);
	}

	if(do_rijk) {
		nrofaux = run_rijk(ref_wfn,pref);
	}
	if(do_rimp2) {
		nrofaux2 = run_rimp2(ref_wfn,pref);
	}

	run_export_mos(ref_wfn,pref);
	run_export_sys(ref_wfn,pref,nrofint,nrofaux,nrofaux2);































/*









    if (dryrun == false)
    {
        double* dumd = new double[aointl+1024];
        int incre = 0;
        Coord   = &(dumd[incre]); incre += 3*nroa;
        charges = &(dumd[incre]); incre +=   nroa;
        mass    = &(dumd[incre]); incre +=   nroa;
        hmat    = &(dumd[incre]); incre += nroao*nroao;
        kmat    = &(dumd[incre]); incre += nroao*nroao;
        smat    = &(dumd[incre]); incre += nroao*nroao;
        Dx      = &(dumd[incre]); incre += nroao*nroao;
        Dy      = &(dumd[incre]); incre += nroao*nroao;
        Dz      = &(dumd[incre]); incre += nroao*nroao;

        MOensA  =  &(dumd[incre]); incre += nroao;
        MOensB  =  &(dumd[incre]); incre += nroao;
        MOsA   = &(dumd[incre]); incre += nroao*nroao;
        MOsB   = &(dumd[incre]); incre += nroao*nroao;


        for (int j=0; j<molecule->natom(); j++) {
            charges[j] = molecule->Z(j);
            mass[j]    = molecule->mass(j);
            for (int i=0; i<3; i++) {
                Coord[3*j+i] = coord->get(j,i);
            }
        }





        if(doTei) {
            count = 0;
            for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
                // Compute quartet
                eri->compute_shell(shellIter);
                // From the quartet get all the integrals
                AOIntegralsIterator intIter = shellIter.integrals_iterator();
                for (intIter.first(); intIter.is_done() == false; intIter.next()) {
                    int p = intIter.i();
                    int q = intIter.j();
                    int r = intIter.k();
                    int s = intIter.l();
                    if( fabs(buffer[intIter.index()]) > cutoff2el) {
                        //	  outfile->Printf("\t(%2d %2d | %2d %2d) = %20.15f\n",
                        //			  p, q, r, s, buffer[intIter.index()]);
                        intval[count] = buffer[intIter.index()];
                        intnums[4*count+0] = p;
                        intnums[4*count+1] = q;
                        intnums[4*count+2] = r;
                        intnums[4*count+3] = s;
 ++count;
                    }
                }
            }
            resort_integrals(intnums, intval,  nrofint, sortcount);
        }
        std::ofstream datf;

        if(doTei) {
            datf.open(sysfn);

            std::cout << "Writing binary data to " <<  "plugin.sys" << "\n";

            double coc[] = {0.,0.,0.};
            double com[] = {0.,0.,0.};
            double tc = 0.;
            double tm = 0.;
            for(int x = 0; x < nroa; x++) {
                tc += charges[x];
                tm += mass[x];
                for(int y = 0; y < 3; y++) {
                    coc[y] += charges[x]*Coord[3*x+y];
                    com[y] += mass[x]*Coord[3*x+y];
                }
            }
            for(int y = 0; y < 3; y++) {
                coc[y] /= tc;
                com[y] /= tm;
            }

            std::cout.precision(12);
            std::cout << "Center of charge is " << coc[0] << " " << coc[1] << " " << coc[2] << "\n";
            std::cout << "Center of mass   is " << com[0] << " " << com[1] << " " << com[2] << "\n";

            //SYSTEM DATA
            datf.write((char *) &nroe, sizeof(int));
            datf.write((char *) &nroao, sizeof(int));
            datf.write((char *) &nroa, sizeof(int));
            datf.write((char *) &nrofint,  sizeof(long long int));
            datf.write((char *) basisName.c_str(), basisName.size());
            datf.write((char *) Coord, sizeof(double)*3*nroa);
            datf.write((char *) charges, sizeof(double)*nroa);
            datf.write((char *) mass,    sizeof(double)*nroa);

            //TWO EL INTEGRAL DATA
            datf.write((char *) sortcount, sizeof(long long int)*4);
            datf.write((char *) intval,    sizeof(double)*nrofint);
            datf.write((char *) intnums,   sizeof(unsigned short)*nrofint*4);

            datf.close();
        }



        delete[] intval;
        delete[] intnums;
    } //end if dyrun





    std::clog<<"\n Starting JKFIT \n";
    auto aux = ref_wfn->get_basisset("JKFIT");

    //calc Metric
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory_J(new IntegralFactory(aux, zero, aux, zero));
    int naux = aux->nbf();
    std::clog<<"naux:"<<naux<<"\n";
    std::shared_ptr<TwoBodyAOInt> Jint (rifactory_J->eri());
    SharedMatrix AOmetric(new Matrix("AO Basis DF Metric", naux, naux));
    double** W = AOmetric->pointer(0);
    const double *Jbuffer = Jint->buffer();
    for (int MU=0; MU < aux->nshell(); ++MU) {
        int nummu = aux->shell(MU).nfunction();
        for (int NU=0; NU <= MU; ++NU) {
            int numnu = aux->shell(NU).nfunction();
            Jint->compute_shell(MU, 0, NU, 0);
            int index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = aux->shell(MU).function_index() + mu;
                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = aux->shell(NU).function_index() + nu;
                    W[omu][onu] = Jbuffer[index];
                    W[onu][omu] = Jbuffer[index];
                }
            }
        }
    }

    SharedMatrix eigvec(new Matrix("eigvecs",naux,naux));
    double** vecp = eigvec->pointer();

    C_DCOPY(naux*(size_t)naux,W[0],1,vecp[0],1);
    double tol = 1.0E-10;
    double* eigval = new double[naux];
    int lwork = naux * 3;
    double* work = new double[lwork];
    int stat = C_DSYEV('v','u',naux,vecp[0],naux,eigval,work,lwork);

    SharedMatrix Jcopy(new Matrix("Jcopy", naux, naux));
    double** Jcopyp = Jcopy->pointer();
    C_DCOPY(naux*(size_t)naux,vecp[0],1,Jcopyp[0],1);

    double max_J = eigval[naux-1];
    int nsig = 0;
    for (int ind=0; ind<naux; ind++) {
        if (eigval[ind] / max_J < tol || eigval[ind] <= 0.0)
            eigval[ind] = 0.0;
        else {
            nsig++;
            eigval[ind] = 1.0 / sqrt(eigval[ind]);
        }
        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(naux, eigval[ind], vecp[ind], 1);
    }
    delete[] eigval;
    C_DGEMM('T','N',naux,naux,naux,1.0,Jcopyp[0],naux,vecp[0],naux,0.0,W[0],naux);

    std::clog<<"\n Metric JKFIT done \n";
    //end metrix


    //calc 3index integrals first round
    //SharedMatrix B(new Matrix("Bso", naux, nbf * nbf));
    //SharedMatrix A(new Matrix("Aso", naux, nbf * nbf));
    //double** Ap = A->pointer();
    //double** Bp = B->pointer();
    //double** Jp = AOmetric->pointer();

    zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(aux,zero,aoBasis,aoBasis));
    std::shared_ptr<TwoBodyAOInt> auxeri(fact->eri());
    buffer = auxeri->buffer();

    count = 0;
    for (int P = 0; P < aux->nshell(); P++) {
        int np = aux->shell(P).nfunction();
        int pstart = aux->shell(P).function_index();
        for (int M = 0; M < aoBasis->nshell(); M++) {
            int nm = aoBasis->shell(M).nfunction();
            int mstart = aoBasis->shell(M).function_index();
            for (int N = 0; N < aoBasis->nshell(); N++) {
                int nn = aoBasis->shell(N).nfunction();
                int nstart = aoBasis->shell(N).function_index();

                auxeri->compute_shell(P,0,M,N);
                for (int p = 0, index = 0; p < np; p++) {
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++, index++) {
                            //Bp[p + pstart][(m + mstart) * nbf + (n + nstart)] = buffer[index];
                            if (buffer[index]  > cutoff2el) {
                                count++;
                                if (count%(long long int)1.25E8 == 0) {
                                    std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::clog<<"JKFIT\nSieve: "<<count*8<<" bytes\n"<<"Total: "<<(long long int)nbf*nbf*naux*8<<" bytes \n\n";
    std::clog<<"\nnbf:"<<nbf<<"\nnaux1:"<<naux;
    //C_DGEMM('N','N',naux, nbf * nbf, naux, 1.0, W[0], naux, Bp[0], nbf * nbf, 0.0,
    //        Ap[0], nbf * nbf);

    // Build numpy and final matrix shape


    std::vector<int> nshape{naux, naux};
    AOmetric->set_numpy_shape(nshape);
 */
	return ref_wfn;


}

}
}  // End Namespaces

inline void swap_ints(long long int a, long long int b, double* intvals, unsigned short* intnums){
	//temp = a
	unsigned short ta = intnums[a*4+0];
	unsigned short tb = intnums[a*4+1];
	unsigned short tc = intnums[a*4+2];
	unsigned short td = intnums[a*4+3];
	double tval = intvals[a];

	//a = b
	intnums[a*4+0] =   intnums[b*4+0];
	intnums[a*4+1] =   intnums[b*4+1];
	intnums[a*4+2] =   intnums[b*4+2];
	intnums[a*4+3] =   intnums[b*4+3];
	intvals[a]     =   intvals[b];

	//b = temp
	intnums[b*4+0] =    ta;
	intnums[b*4+1] =    tb;
	intnums[b*4+2] =    tc;
	intnums[b*4+3] =    td;
	intvals[b]     =    tval;

}

void resort_integrals(unsigned short int*  intnums, double* intval,  long long int nrofint, long long int* sortcount ){
	//RESORT INTEGRALS

	//STEP 1: BRING to basis types
	for(long long int x = 0; x < nrofint; x++) {
		unsigned short a,b,c,d;

		//INPUT AUCH  CHEMIKER ??
		a      = intnums[x*4+0];//a
		b      = intnums[x*4+1];//b
		c      = intnums[x*4+2];//c
		d      = intnums[x*4+3];//d

		//TYP IIb: reorder to IIa
		if(a==b && b==d && c!=d) {
			intnums[x*4+0] = a; //a
			intnums[x*4+1] = a; //b
			intnums[x*4+2] = a; //c
			intnums[x*4+3] = c; //d
		}

		//TYP IIc: reorder to IIa
		if(a!=b && a==c && c==d) {
			intnums[x*4+0] = a; //a
			intnums[x*4+1] = a; //b
			intnums[x*4+2] = a; //c
			intnums[x*4+3] = b; //d
		}

		//TYP IId: reorder to IIa
		if(a!=b && b==c && c==d) {
			intnums[x*4+0] = b; //a
			intnums[x*4+1] = b; //b
			intnums[x*4+2] = b; //c
			intnums[x*4+3] = a; //d
		}

		//TYP IIIc: reorder to IIIb
		if(a==d && b==c && a!=b) {
			intnums[x*4+0] = a; //a
			intnums[x*4+1] = b; //b
			intnums[x*4+2] = d; //c
			intnums[x*4+3] = c; //d
		}


		//TYP IVc:  reorder IVb
		if(b==c && b!=a && b!=d && a!=d) {
			intnums[x*4+0] = b; //a
			intnums[x*4+1] = a; //b
			intnums[x*4+2] = c; //c
			intnums[x*4+3] = d; //d
		}

		//TYP IVd:  reorder to IVa
		if(c==d && c!=a && c!=b && a!=b) {
			intnums[x*4+0] = c; //a
			intnums[x*4+1] = c; //b
			intnums[x*4+2] = a; //c
			intnums[x*4+3] = b; //d
		}

		//TYP IVe:  reorder IVb
		if(a==d && a!=b && a!=c && b!=c) {
			intnums[x*4+0] = a; //a
			intnums[x*4+1] = b; //b
			intnums[x*4+2] = d; //c
			intnums[x*4+3] = c; //d
		}
		//TYP IVf: perm_all,  reorder IVb
		if(b==d && b!=a && b!=c && a!=c) {
			intnums[x*4+0] = b; //a
			intnums[x*4+1] = a; //b
			intnums[x*4+2] = d; //c
			intnums[x*4+3] = c; //d
		}
	}

	//STEP2: BRING PERM1 to front
	long long int sorted = 0;
	for(long long int x = 0; x < nrofint; x++) {
		unsigned short a,b,c,d;

		//INPUT AUCH  CHEMIKER ??
		a      = intnums[x*4+0];//a
		b      = intnums[x*4+1];//b
		c      = intnums[x*4+2];//c
		d      = intnums[x*4+3];//d
		if(a==b && b==c && c==d) {
			swap_ints(x, sorted,  intval,  intnums);
			sorted++;
		}
	}
	std::clog << sorted << " integrals after search for perm_1\n";
	sortcount[0] = sorted;
	//STEP3: BRING PERM1_5 thereafter

	for(long long int x = sorted; x < nrofint; x++) {
		unsigned short a,b,c,d;

		//INPUT AUCH  CHEMIKER ??
		a      = intnums[x*4+0];//a
		b      = intnums[x*4+1];//b
		c      = intnums[x*4+2];//c
		d      = intnums[x*4+3];//d
		if(a==b && c==d) {
			swap_ints(x, sorted,  intval,  intnums);
			sorted++;
		}
	}
	std::clog << sorted << " integrals after search for perm_15\n";
	sortcount[1] = sorted;


	//STEP4: BRING PERM1234
	for(int x = sorted; x < nrofint; x++) {
		int a,b,c,d;

		//INPUT AUCH  CHEMIKER ??
		a      = intnums[x*4+0];//a
		b      = intnums[x*4+1];//b
		c      = intnums[x*4+2];//c
		d      = intnums[x*4+3];//d
		if(a==c && b==d) {
			swap_ints(x, sorted,  intval,  intnums);
			sorted++;
		}
	}
	std::clog << sorted << " integrals after search for perm_1234\n";
	sortcount[2] = sorted;

	//STEP5: BRING PERM1256
	for(int x = sorted; x < nrofint; x++) {
		int a,b,c,d;

		//INPUT AUCH  CHEMIKER ??
		a      = intnums[x*4+0];//a
		b      = intnums[x*4+1];//b
		c      = intnums[x*4+2];//c
		d      = intnums[x*4+3];//d
		if(a==b) {
			swap_ints(x, sorted,  intval,  intnums);
			sorted++;
		}
	}
	std::clog << sorted << " integrals after search for perm_1256\n";
	sortcount[3] = sorted;
	std::clog << nrofint << " integrals in total\n";
	std::clog << "-------------------------------------------------------------------------------\n";

}
