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

void resort_integrals(unsigned short int*  intnums, double* intval,  long long int nrofint, long long int* sortcount );


namespace psi {

namespace syshfwplug  {

extern "C"
int read_options(std::string name, Options &options)
{
	if (name == "SYSHFWPLUG"|| options.read_globals()) {
		/*- The amount of information printed
		    to the output file -*/
		options.add_int("PRINT", 1);
		/*- Whether to compute two-electron integrals -*/
		options.add_bool("DO_TEI", true);
		options.add_bool("DRYRUN", false);
		options.add_bool("BUILD_RIJK ", false);
		options.add_bool("BUILD_RIMP2", false);
		options.add_str("PREF","vergessen-wa");

	}

	return true;
}

extern "C"
SharedMatrix syshfwplug(SharedWavefunction ref_wfn, Options &options)
{
	// Grab options from the options object
	int print = options.get_int("PRINT");
	int doTei = options.get_bool("DO_TEI");
	std::string pref = options.get_str("PREF");
	bool dryrun = options.get_bool("DRYRUN");

	std::cout << "Prefix " << pref << "\n";
	std::cout << "Write sys " << doTei << "\n";

	std::string sysfn = pref + ".sys";
	std::string hfawfn = pref + ".ahfw";
	std::string hfbwfn = pref + ".bhfw";



	std::shared_ptr<Molecule> molecule = ref_wfn->molecule();

	// Form basis object:
	std::shared_ptr<BasisSet> aoBasis = ref_wfn->basisset();

	// The integral factory oversees the creation of integral objects
	std::shared_ptr<IntegralFactory> integral(new IntegralFactory(aoBasis, aoBasis, aoBasis, aoBasis));

	// N.B. This should be called after the basis has been built, because the geometry has not been
	// fully initialized until this time.

	//shared_ptr<Vector3> bal = molecule->xyz(0); ;
	//    bal = molecule->xyz(0);


	const double cutoff2el = 1.e-12;

	molecule->print();
	std::cout << "Nr of atoms " <<  molecule->natom() << "\n";
	std::shared_ptr<Matrix> coord (new Matrix(molecule->geometry()));


	std::cout << "DEBUG";
	// Build a dimension object with a single dim since this is AO's

	int nbf = aoBasis->nbf();

	double nucrep = molecule->nuclear_repulsion_energy();
	outfile->Printf("\n    Nuclear repulsion energy: %16.8f\n\n", nucrep);

	// The matrix factory can create matrices of the correct dimensions...
	std::shared_ptr<MatrixFactory> factory(new MatrixFactory);
	factory->init_with(1, &nbf, &nbf);

	// Form the one-electron integral objects from the integral factory
	std::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
	std::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
	std::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
	std::shared_ptr<OneBodyAOInt> ecpOBI(integral->ao_ecp());
	std::shared_ptr<OneBodyAOInt> dip(integral->ao_dipole());

	// Form the one-electron integral matrices from the matrix factory
	std::shared_ptr<Matrix> sMat(factory->create_matrix("Overlap"));
	std::shared_ptr<Matrix> tMat(factory->create_matrix("Kinetic"));
	std::shared_ptr<Matrix> vMat(factory->create_matrix("Potential"));
	std::shared_ptr<Matrix> hMat(factory->create_matrix("One Electron Ints"));
	std::shared_ptr<Matrix> ecpMat(factory->create_matrix("EPC Matrix "));



	std::vector<SharedMatrix> dipole;
	dipole.push_back(SharedMatrix(new Matrix("AO Mux", aoBasis->nbf(), aoBasis->nbf())));
	dipole.push_back(SharedMatrix(new Matrix("AO Muy", aoBasis->nbf(), aoBasis->nbf())));
	dipole.push_back(SharedMatrix(new Matrix("AO Muz", aoBasis->nbf(), aoBasis->nbf())));

	dip->compute(dipole);

	// Compute the one electron integrals, telling each object where to store the result


	// Form h = T + V by first cloning T and then adding V

	hMat->copy(tMat);
	hMat->add(vMat);
	hMat->add(ecpMat);


	sOBI->compute(sMat);
	std::shared_ptr<Matrix> Calpha;
	std::shared_ptr<Matrix> Cbeta;
	//Coeff.
	if (ref_wfn->reference_energy() != 0.0)
	{
		Calpha = ref_wfn->Ca();   //last index over MOS!!!!!!
		Cbeta  = ref_wfn->Cb();
		hMat   = ref_wfn->H();
	}
	else
	{
		Calpha = std::make_shared<Matrix>(aoBasis->nbf(),aoBasis->nbf());
		Cbeta  = std::make_shared<Matrix>(aoBasis->nbf(),aoBasis->nbf());

		if (dryrun==false) {
			tOBI->compute(tMat);
			vOBI->compute(vMat);
			ecpOBI->compute(ecpMat);
			hMat->copy(tMat);
			hMat->add(vMat);
			hMat->add(ecpMat);
		}

	}


	long long int count=0;

	// Now, the two-electron integrals
	std::shared_ptr<TwoBodyAOInt> eri(integral->eri());
	// The buffer will hold the integrals for each shell, as they're computed
	const double *buffer = eri->buffer();
	// The iterator conveniently lets us iterate over functions within shells
	AOShellCombinationsIterator shellIter = integral->shells_iterator();
	if(doTei) {
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
				}
			}
		}
	}

	long long int nrofint = count;       //nr of nonzero to electron integrals
	long long int sortcount[4];  //boundaries for different permutation patterns;
	int nroao = aoBasis->nbf();
	int nroa  = molecule->natom();
	int nroe  = ref_wfn->nalpha() + ref_wfn->nbeta();

	double*              intval;        //two electron integrals
	unsigned short int*  intnums;       //two electron indecies
	double*              Coord;         //atomic ccordinates
	double*              charges;       //atomic charges  (all set to the same value so far)
	double*              mass;          //atomic masses
	double*              hmat;          //one electron hamiltonian
	double*              kmat;          //kinetic energy
	double*              smat;          //overlapp
	double*              Dx;            //dipole x
	double*              Dy;            //dipole y
	double*              Dz;            //dipole z

	double* MOensA;   //nroao
	double* MOensB;   //nroao
	double* MOsA;     //nroao*nroao
	double* MOsB;     //nroao*nroao

	int aointl = nroao*nroao*8+nroa*5+2*nroao;
	std::cout << "Need " <<  aointl*8 <<   " Bytes  for one el ints, coords, and charges and MOs.\n";
	std::cout << "Need " <<  nrofint*8 <<  " Bytes  for two el indices and the same for two el values.\n";
	std::cout<<"\nORBBAS\nSieve: "<<nrofint*8<<" bytes\n"<<"Total: "<<(long long int) nbf*nbf*nbf*nbf*8<<" bytes \n\n";


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

		intval = new double[nrofint];
		intnums = new unsigned short[nrofint*4+4];
		for (int j=0; j<molecule->natom(); j++) {
			charges[j] = molecule->Z(j);
			mass[j]    = molecule->mass(j);
			for (int i=0; i<3; i++) {
				Coord[3*j+i] = coord->get(j,i);
			}
		}

		for(int x = 0; x < nroao; x++) {
			for(int y = 0; y < nroao; y++) {
				hmat[x*nroao+y] = hMat->get(x,y);
				kmat[x*nroao+y] = tMat->get(x,y);
				smat[x*nroao+y] = sMat->get(x,y);
				Dx[x*nroao+y] = -dipole[0]->get(x,y);
				Dy[x*nroao+y] = -dipole[1]->get(x,y);
				Dz[x*nroao+y] = -dipole[2]->get(x,y);
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

			std::clog << "Writing binary data to " <<  "plugin.sys" << "\n";

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

			std::clog.precision(12);
			std::clog << "Center of charge is " << coc[0] << " " << coc[1] << " " << coc[2] << "\n";
			std::clog << "Center of mass   is " << com[0] << " " << com[1] << " " << com[2] << "\n";

			//SYSTEM DATA
			datf.write((char *) &nroe, sizeof(int));
			datf.write((char *) &nroao, sizeof(int));
			datf.write((char *) &nroa, sizeof(int));
			datf.write((char *) &nrofint,  sizeof(long long int));
			datf.write((char *) Coord, sizeof(double)*3*nroa);
			datf.write((char *) charges, sizeof(double)*nroa);
			datf.write((char *) mass,    sizeof(double)*nroa);

			//ONEL EL INTEGRAL DATA
			datf.write((char * ) hmat, sizeof(double)*nroao*nroao);
			datf.write((char * ) kmat, sizeof(double)*nroao*nroao);
			datf.write((char * ) smat, sizeof(double)*nroao*nroao);
			datf.write((char * ) Dx, sizeof(double)*nroao*nroao);
			datf.write((char * ) Dy, sizeof(double)*nroao*nroao);
			datf.write((char * ) Dz, sizeof(double)*nroao*nroao);

			//TWO EL INTEGRAL DATA
			datf.write((char *) sortcount, sizeof(long long int)*4);
			datf.write((char *) intval,    sizeof(double)*nrofint);
			datf.write((char *) intnums,   sizeof(unsigned short)*nrofint*4);

			datf.close();
		}

		for(int x = 0; x < nroao; x++)
			MOensA[x] = MOensB[x] = 0.;

		for(int x = 0; x < nroao; x++) {
			for(int y = 0; y < nroao; y++) {
				MOsA[x*nroao+y] = Calpha->get(y,x);
				MOsB[x*nroao+y] = Cbeta->get(y,x);
			}
		}

		std::clog << "A " << MOsA[0] << " B " << MOsB[0] << "\n";
		std::clog << "Writing MOs to " << hfawfn << " & " << hfbwfn << "\n";
		datf.open(hfawfn);
		datf.write((char *)  &nroao, sizeof(int));
		datf.write((char * ) MOensA, sizeof(double)*nroao);
		datf.write((char * ) MOsA, sizeof(double)*nroao*nroao);
		datf.close();

		datf.open(hfbwfn);
		datf.write((char *)  &nroao, sizeof(int));
		datf.write((char * ) MOensB, sizeof(double)*nroao);
		datf.write((char * ) MOsB, sizeof(double)*nroao*nroao);
		datf.close();

		delete[] intval;
		delete[] intnums;
	} //end if dyrun





  std::cout<<"\n Starting JKFIT \n";
	auto aux = ref_wfn->get_basisset("JKFIT");

	//calc Metric
	std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
	std::shared_ptr<IntegralFactory> rifactory_J(new IntegralFactory(aux, zero, aux, zero));
	int naux = aux->nbf();
	std::cout<<"naux:"<<naux<<"\n";
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

	std::cout<<"\n Metric JKFIT done \n";
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
							if (buffer[index]  > cutoff2el){
								count++;
								if (count%(long long int)1.25E8 == 0){
									std::cout<<count/(long long int)1.25E8<<'\t'<<std::flush;
								}
							}
						}
					}
				}
			}
		}
	}
	std::cout<<"JKFIT\nSieve: "<<count*8<<" bytes\n"<<"Total: "<<(long long int)nbf*nbf*naux*8<<" bytes \n\n";
	std::cout<<"\nnbf:"<<nbf<<"\nnaux1:"<<naux;
	//C_DGEMM('N','N',naux, nbf * nbf, naux, 1.0, W[0], naux, Bp[0], nbf * nbf, 0.0,
	//        Ap[0], nbf * nbf);

	// Build numpy and final matrix shape


	std::vector<int> nshape{naux, naux};
	AOmetric->set_numpy_shape(nshape);
	return AOmetric;


	//auto aux2 = ref_wfn->get_basisset("RIFIT");
	//std::clog<<aux2->nbf();



	//stolen from CDLK.cc

	/*
	   auto sieve = std::make_shared<ERISieve>(aoBasis, 0.0);
	   long int ncholesky_;

	   double cholesky_tolerance_ = 0.001;
	   size_t memory_ = 500;

	   std::clog<<"1";


	   auto integral2 = std::make_shared<IntegralFactory>(aoBasis,aoBasis,aoBasis,aoBasis);

	   auto Ch = std::make_shared<CholeskyERI>(std::shared_ptr<TwoBodyAOInt>(integral2->eri()),0.0,cholesky_tolerance_,memory_);
	   std::clog<<"1.2";


	   int ntri = sieve->function_pairs().size();
	   std::clog<<"2";

	   Ch->choleskify();
	   std::clog<<"3";


	   ncholesky_  = Ch->Q();
	   size_t three_memory = ncholesky_ * ntri;
	   int nbf_int = aoBasis->nbf();
	   std::clog<<"2";

	   std::shared_ptr<Matrix> L = Ch->L();
	   double ** Lp = L->pointer();

	   auto Qmn_ = std::make_shared<Matrix>("Qmn (CD Integrals)", ncholesky_ , ntri);
	   double** Qmnp = Qmn_->pointer();

	   const std::vector<long int>& schwarz_fun_pairs = sieve->function_pairs_reverse();

	   for (size_t mu = 0; mu < nbf_int; mu++) {
	    for (size_t nu = mu; nu < nbf_int; nu++) {
	        if ( schwarz_fun_pairs[nu*(nu+1)/2+mu] < 0 ) continue;
	        for (long int P = 0; P < ncholesky_; P++) {
	            Qmnp[P][schwarz_fun_pairs[nu*(nu+1)/2+mu]] = Lp[P][mu * nbf_int + nu];
	        }
	    }
	   }


	 */


















	// Obtain the Wavefunction from globals that we set python-side
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
