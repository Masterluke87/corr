#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "ops_rhf.h"
#include <omp.h>
#include <libint2.hpp>


//Extern Functions
extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);
extern void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
extern void transform_MOs(int nroao, double *MOs, double* tmat, double* tmpvec);


void calculate_libint_oei(std::vector<libint2::Atom> &atoms,libint2::BasisSet &obs,double* zeff,
                          double* Hmat,double* Tmat,double* Smat,double* Vmat,
                          double* Hmat_libint,double* Tmat_libint,double* Smat_libint,double* Vmat_libint,
                          double* Hmat_trans,double* Tmat_trans,double* Smat_trans,double* Vmat_trans)
{



	libint2::Engine s_engine(libint2::Operator::overlap,obs.max_nprim(),obs.max_l());
	libint2::Engine t_engine(libint2::Operator::kinetic,obs.max_nprim(),obs.max_l());
	libint2::Engine v_engine(libint2::Operator::nuclear,obs.max_nprim(),obs.max_l());

	std::vector<std::pair<double, std::array<double, 3> > > q(atoms.size());
	for (size_t i = 0; i < atoms.size(); i++) {
		q.push_back({zeff[i],{{atoms[i].x,atoms[i].y,atoms[i].z}}});
	}
	v_engine.set_params(q);

	int nroao = obs.nbf();

	auto shell2bf = obs.shell2bf();
	std::vector<int> one_shift;
	std::vector<int> two_shift;
	int one_size = 0;
	int two_size = 0;
	const auto& buf_vec_t = t_engine.results();
	const auto& buf_vec_s = s_engine.results();
	const auto& buf_vec_v = v_engine.results();

	for(auto s1=0; s1!=obs.size(); ++s1) {
		for(auto s2=0; s2!=obs.size(); ++s2) {
			t_engine.compute(obs[s1], obs[s2]);
			s_engine.compute(obs[s1], obs[s2]);
			v_engine.compute(obs[s1], obs[s2]);

			auto ints_shellset_t = buf_vec_t[0];
			auto ints_shellset_s = buf_vec_s[0];
			auto ints_shellset_v = buf_vec_v[0];

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
			if (obs[s1].contr[0].l == 5) {
				one_shift = {+10,+7,+4,+1,-2,-5,-5,-4,-3,-2,-1};
				one_size  = 11;
			}
			if (obs[s2].contr[0].l == 5) {
				two_shift = {+10,+7,+4,+1,-2,-5,-5,-4,-3,-2,-1};
				two_size  = 11;
			}
			// integrals are packed into ints_shellset in row-major (C) form
			// this iterates over integrals in this order
			for(auto f1=0; f1!=n1; ++f1)
				for(auto f2=0; f2!=n2; ++f2) {
					Tmat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_t[f1*n2+f2];
					Smat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_s[f1*n2+f2];
					Vmat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_v[f1*n2+f2];

					Tmat_trans[(bf1+f1)*obs.nbf() + (bf2+f2)] = Tmat[(bf1+f1+one_shift[f1])*obs.nbf() + (bf2+f2+two_shift[f2])];
					Smat_trans[(bf1+f1)*obs.nbf() + (bf2+f2)] = Smat[(bf1+f1+one_shift[f1])*obs.nbf() + (bf2+f2+two_shift[f2])];
					Vmat_trans[(bf1+f1)*obs.nbf() + (bf2+f2)] = Vmat[(bf1+f1+one_shift[f1])*obs.nbf() + (bf2+f2+two_shift[f2])];

					Hmat_trans[(bf1+f1)*obs.nbf() + (bf2+f2)] = Hmat[(bf1+f1+one_shift[f1])*obs.nbf() + (bf2+f2+two_shift[f2])];
				}
		}
	}
	for (size_t i = 0; i < nroao*nroao; i++) {
		Hmat_libint[i] = Tmat_libint[i] + Vmat_libint[i];
	}

	double tmatdiff = 0;
	double smatdiff = 0;
	double vmatdiff = 0;
	double hmatdiff = 0;

	double tmatmax = 0;
	double smatmax = 0;
	double vmatmax = 0;
	double hmatmax = 0;


	for (size_t i = 0; i < nroao*nroao; i++) {
		tmatdiff += fabs(Tmat_trans[i] - Tmat_libint[i]);
		smatdiff += fabs(Smat_trans[i] - Smat_libint[i]);
		vmatdiff += fabs(Vmat_trans[i] - Vmat_libint[i]);
		hmatdiff += fabs(Hmat_trans[i] - Hmat_libint[i]);

		if (fabs(Tmat_trans[i] - Tmat_libint[i])>tmatmax)
			tmatmax = fabs(Tmat_trans[i] - Tmat_libint[i]);
		if (fabs(Smat_trans[i] - Smat_libint[i]) > smatmax)
			smatmax = fabs(Smat_trans[i] - Smat_libint[i]);
		if (fabs(Vmat_trans[i] - Vmat_libint[i]) > vmatmax)
			vmatmax = fabs(Vmat_trans[i] - Vmat_libint[i]);
		if (fabs(Hmat_trans[i] - Hmat_libint[i]) > hmatmax)
			hmatmax = fabs(Hmat_trans[i] - Hmat_libint[i]);


	}

	std::cout << "Capability check (did the transformation work?): " << '\n';
	std::cout << "tmatdiff:" <<tmatdiff << '\n';
	std::cout << "smatdiff:" <<smatdiff << '\n';
	std::cout << "vmatdiff:" <<vmatdiff << '\n';
	std::cout << "hmatdiff:" <<hmatdiff << '\n';

	std::cout << "tmatmax :" <<tmatmax << '\n';
	std::cout << "smatmax :" <<smatmax << '\n';
	std::cout << "vmatmax :" <<vmatmax << '\n';
	std::cout << "hmatmax :" <<hmatmax << '\n';

	if (smatmax < 1E-10 && tmatmax<1E-10 && vmatmax<1E-10) {
		std::cout << "It seems the transformation worked" << '\n';
	}




}










/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_r_ab                                                       */
/*                                                                               */
/* caluculates the ion core coulomb repulsion                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_r_ab(int a, int b, double* coord){
	double r = sqrt(pow(coord[3*a+0]-coord[3*b+0],2)+
	                pow(coord[3*a+1]-coord[3*b+1],2)+
	                pow(coord[3*a+2]-coord[3*b+2],2));
	return(r);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_ion_rep                                                    */
/*                                                                               */
/* calculates the ion core coulomb repulsion                                     */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_ion_rep(int nroa, double* coord, double* charges){
	double ion_rep = 0.;
	for(int x = 0; x < nroa; x++) {
		for(int y = x+1; y < nroa; y++) {
			ion_rep += charges[x]*charges[y]/calc_r_ab(x,y,coord);
		}
	}

	return(ion_rep);
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_center_of_mass                                             */
/*                                                                               */
/* calculates the ion core center of mass                                         */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass){
	for(int x = 0; x < 3; x++) center_of_mass[x] = 0.;

	double mass_tot = 0.;

	for(int x = 0; x < nroa; x++) {
		mass_tot += mass[x];
		for(int i = 0; i < 3; i++)
			center_of_mass[i] += mass[x] * coord[3*x+i];
	}

	for(int i = 0; i < 3; i++)
		center_of_mass[i] /= mass_tot;

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_mu_core                                                    */
/*                                                                               */
/* calculates the ion core dipole moment  at point                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_mu_core(int nroa, double* coord, double* charges, double* point,
                  double* mu_core){
	for(int x = 0; x < 3; x++) mu_core[x] = 0.;


	for(int x = 0; x < nroa; x++) {
		for(int i = 0; i < 3; i++)
			mu_core[i] += charges[x] * (coord[3*x+i] - point[i]);
	}

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_S12                                                        */
/*                                                                               */
/* calculate S^-1/2 for symmetric orthogonalisation                              */
/* returns mimum eigenvalue of S                                                 */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_S12(int nroao, double* Smat, double* Som12){

	double* tmpspace = new double[3*nroao*nroao];
	int inc=0;
	double* tmpmat = &(tmpspace[inc]); inc+=nroao*nroao;
	double* tmpvecs = &(tmpspace[inc]); inc+=nroao*nroao;
	double* tmpvals = &(tmpspace[inc]);

	diag_mat(nroao, Smat, tmpvals, tmpvecs);

	double min_val =10.;
	for(int x = 0; x < nroao; x++) {
		if(min_val > tmpvals[x]) min_val = tmpvals[x];
	}

	for(int x = 0; x <  nroao*nroao; x++) {
		Som12[x] = 0.; tmpmat[x] = 0.;
	}
	for(int x = 0; x <  nroao; x++)
		Som12[x*nroao+x] = 1./sqrt(tmpvals[x]);
	for(int x = 0; x  < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			for(int z = 0; z < nroao; z++) //Adjungiert !!!!!!
				tmpmat[x*nroao+y] +=  Som12[z*nroao+y] * tmpvecs[z*nroao+x];
		}
	}
	for(int x = 0; x  < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			Som12[x*nroao+y] = 0.;
			for(int z = 0; z < nroao; z++)
				Som12[x*nroao+y] += tmpvecs[z*nroao+y] * tmpmat[x*nroao+z];
		}
	}
	delete[] tmpspace;
	return(min_val);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_Sp12                                                       */
/*                                                                               */
/* calculate S^+1/2 for symmetric orthogonalisation of P                         */
/* returns mimum eigenvalue of S                                                 */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_Sp12(int nroao, double* Smat, double* Sop12, double* tmpmat,
                 double* tmpvecs, double* tmpvals){

	diag_mat(nroao, Smat, tmpvals, tmpvecs);

	double min_val =10.;
	for(int x = 0; x < nroao; x++) {
		if(min_val > tmpvals[x]) min_val = tmpvals[x];
	}

	for(int x = 0; x <  nroao*nroao; x++) {
		Sop12[x] = 0.; tmpmat[x] = 0.;
	}
	for(int x = 0; x <  nroao; x++)
		Sop12[x*nroao+x] = sqrt(tmpvals[x]);
	for(int x = 0; x  < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			for(int z = 0; z < nroao; z++) //Adjungiert !!!!!!
				tmpmat[x*nroao+y] +=  Sop12[z*nroao+y] * tmpvecs[z*nroao+x];
		}
	}
	for(int x = 0; x  < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			Sop12[x*nroao+y] = 0.;
			for(int z = 0; z < nroao; z++)
				Sop12[x*nroao+y] += tmpvecs[z*nroao+y] * tmpmat[x*nroao+z];
		}
	}

	return(min_val);
}

/*******************************************************************************
 * This routine diagonalises the fock operator                                 *
 *                                                                             *
 ******************************************************************************/

void diag_Fmat(int nroao, double* Fmat, double* MOs,
               double* MOens, double* Som12, double* tmpmat){
	symmortho_mat(nroao,Fmat,Som12,tmpmat);
	diag_mat(nroao,Fmat,MOens,MOs);
	transform_MOs(nroao, MOs, Som12,tmpmat);
}


/*******************************************************************************
 * This routine builds the new Pmat in a damped SCF                            *
 *                                                                             *
 ******************************************************************************/

double build_Pmat_dscf(int nroao, int nroe, double* Pmat, double* Pmat_old,
                       double* MOs, double damp){
	double max_diff = 0.;

	for(int x = 0; x < nroao*nroao; x++) Pmat[x] = 0.;

	for(int e = 0; e < nroe/2; e++) {
		for(int x = 0; x < nroao; x++) {
			for(int y = 0; y < nroao; y++) {
				Pmat[x*nroao+y] += 2.*MOs[e*nroao+y]*MOs[e*nroao+x];
			}
		}
	}

	for(int x = 0; x < nroao*nroao; x++) {
		if(fabs(Pmat_old[x] - Pmat[x]) > max_diff) max_diff = fabs(Pmat_old[x] - Pmat[x]);
		Pmat[x] = (1.-damp)*Pmat[x] + damp*Pmat_old[x];
		Pmat_old[x] = Pmat[x];
	}

	return(max_diff);
}

/*******************************************************************************
 * Permutations for fock matrix build                                          *
 *                                                                             *
 ******************************************************************************/


void perm_all(unsigned short a, unsigned short b,
              unsigned short c, unsigned short d, double integral,
              double* Fmat, double* Pmat, int nroao){
	//J-Terme
	Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
	Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
	Fmat[b*nroao+a] += Pmat[c*nroao+d]*integral; //3
	Fmat[b*nroao+a] += Pmat[d*nroao+c]*integral; //4
	Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5
	Fmat[d*nroao+c] += Pmat[a*nroao+b]*integral; //6
	Fmat[c*nroao+d] += Pmat[b*nroao+a]*integral; //7
	Fmat[d*nroao+c] += Pmat[b*nroao+a]*integral; //8

	//K-Terme
	Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
	Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
	Fmat[b*nroao+c] -= 0.5*Pmat[a*nroao+d]*integral; //3
	Fmat[b*nroao+d] -= 0.5*Pmat[a*nroao+c]*integral; //4
	Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
	Fmat[d*nroao+a] -= 0.5*Pmat[c*nroao+b]*integral; //6
	Fmat[c*nroao+b] -= 0.5*Pmat[d*nroao+a]*integral; //7
	Fmat[d*nroao+b] -= 0.5*Pmat[c*nroao+a]*integral; //8
}

void perm_1234(unsigned short a, unsigned short b,
               unsigned short c, unsigned short d, double integral,
               double* Fmat, double* Pmat, int nroao){
	//J-Terme
	Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
	Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
	Fmat[b*nroao+a] += Pmat[c*nroao+d]*integral; //3
	Fmat[b*nroao+a] += Pmat[d*nroao+c]*integral; //4

	//K-Terme
	Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
	Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
	Fmat[b*nroao+c] -= 0.5*Pmat[a*nroao+d]*integral; //3
	Fmat[b*nroao+d] -= 0.5*Pmat[a*nroao+c]*integral; //4
}


void perm_1256(unsigned short a, unsigned short b,
               unsigned short c, unsigned short d, double integral,
               double* Fmat, double* Pmat, int nroao){
	//J-Terme
	Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
	Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
	Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5
	Fmat[d*nroao+c] += Pmat[a*nroao+b]*integral; //6

	//K-Terme
	Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
	Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
	Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
	Fmat[d*nroao+a] -= 0.5*Pmat[c*nroao+b]*integral; //6
}

void perm_15(unsigned short a, unsigned short b,
             unsigned short c, unsigned short d, double integral,
             double* Fmat, double* Pmat, int nroao){
	//J-Terme
	Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
	Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5

	//K-Terme
	Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
	Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
}

void perm_1(unsigned short a, unsigned short b,
            unsigned short c, unsigned short d, double integral,
            double* Fmat, double* Pmat, int nroao){
	//J-Terme
	Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
	//K-Terme
	Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
}

/*******************************************************************************
 * Fock matrix build                                                           *
 *                                                                             *
 ******************************************************************************/
void build_Fmat(int nroao, double* Fmat, double* Pmat, double* Hmat,double* intval, unsigned short* intnums, long long int* sortcount, long long int nrofint){
	for(int mu = 0; mu < nroao; mu++) {
		for(int nu = 0; nu < nroao; nu++) {
			Fmat[mu*nroao+nu] = Hmat[mu*nroao+nu];
		}
	}

	//cannot get to work a "global reduction" -> after ever

	//PERM 1
	#pragma omp parallel for reduction(+:Fmat[:nroao*nroao])
	for(long long int x = 0; x < sortcount[0]; x++)
		perm_1(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);

	//PERM_15
	#pragma omp parallel for reduction(+:Fmat[:nroao*nroao])
	for(long long int x = sortcount[0]; x < sortcount[1]; x++) {
		perm_15(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);
	}

	//PERM_1234
	#pragma omp parallel for reduction(+:Fmat[:nroao*nroao])
	for(long long int x = sortcount[1]; x < sortcount[2]; x++)
		perm_1234(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);


	//PERM_1256
	#pragma omp parallel for reduction(+:Fmat[:nroao*nroao])
	for(long long int x = sortcount[2]; x < sortcount[3]; x++)
		perm_1256(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);


	//PERM_ALL
	#pragma omp parallel for reduction(+:Fmat[:nroao*nroao])
	for(long long int x = sortcount[3]; x < nrofint; x++)
		perm_all(intnums[x*4+0],intnums[x*4+1],intnums[x*4+2],intnums[x*4+3],intval[x],Fmat,Pmat,nroao);
}


/*******************************************************************************
 *  Calculation of electronic HF-energy                                        *
 *                                                                             *
 ******************************************************************************/

double Calc_e_el(int nroao, double* Fmat, double* Pmat, double* Hmat){
	double E_tot_el = 0.;
	for(int x = 0; x < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			E_tot_el += Pmat[x*nroao+y]*(Fmat[y*nroao+x] + Hmat[y*nroao+x]);
		}
	}
	return(E_tot_el/2);
}
/*******************************************************************************
 * calculation of one electrion expc. values                                   *
 *                                                                             *
 ******************************************************************************/

double   calc_op_1el(int nroao, double* opmat, double* Pmat){
	double op = 0.;
	for(int x = 0; x < nroao; x++) {
		for(int y = 0; y < nroao; y++) {
			op += Pmat[x*nroao+y]*opmat[y*nroao+x];
		}
	}
	return(op);
}


void calculate_libint_tei(std::vector<libint2::Atom> &atoms,libint2::BasisSet &obs,
                          long long int &nrofint,double** intval,unsigned short** intnums,long long int * sortcount)
{
	libint2::Engine eri_engine(libint2::Operator::coulomb,obs.max_nprim(),obs.max_l());
	auto shell2bf = obs.shell2bf();
	const auto& buf_vec_eri = eri_engine.results();
	long long int count = 0;
	double integral = 0.0;

	double* Schwarz = new double[obs.size()*obs.size()];
	double max_in_shell = 0.0;

	for(auto s2=0; s2!=obs.size(); ++s2)
		for(auto s1=s2; s1!=obs.size(); ++s1)
			for(auto s4=0; s4!=obs.size(); ++s4)
				for(auto s3=s4; s3!=obs.size(); ++s3)  {
					eri_engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
					auto ints_shellset_eri = buf_vec_eri[0];

					if (ints_shellset_eri!=NULL) {
						max_in_shell = 0.0;

						auto bf1 = shell2bf[s1]; // first basis function in first shell
						auto n1 = obs[s1].size(); // number of basis functions in first shell
						auto bf2 = shell2bf[s2]; // first basis function in second shell
						auto n2 = obs[s2].size(); // number of basis functions in second shell
						auto bf3 = shell2bf[s3]; // first basis function in first shell
						auto n3 = obs[s3].size(); // number of basis functions in first shell
						auto bf4 = shell2bf[s4]; // first basis function in second shell
						auto n4 = obs[s4].size(); // number of basis functions in second shell

						for(auto f1=0; f1!=n1; ++f1)
							for(auto f2=0; f2!=n2; ++f2)
								for(auto f3=0; f3!=n3; ++f3)
									for(auto f4=0; f4!=n4; ++f4) {
										//Tmat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_t[f1*n2+f2]
										integral = ints_shellset_eri[f1*n2*n3*n4+f2*n3*n4+f3*n4+f4];
										if (std::fabs(integral)>max_in_shell) {
											max_in_shell=sqrt(std::fabs(integral));
										}

										if (((bf1+f1) >= (bf2+f2)) && ((bf3+f3) >= (bf4+f4)) && (((bf1+f1)*(bf1+f1+1)/2 +(bf2+f2))>=((bf3+f3)*(bf3+f3+1)/2+(bf4+f4))))
										{
											if (std::fabs(integral)>1E-12)
												count++;
										}
									}
						if (s1==s3 && s2 == s4) {
							Schwarz[s1*obs.size()+s2] = max_in_shell;
						}
					}
				}

	*intval        = new double[count];                       //two electron integrals
	*intnums       = new unsigned short[count*4];             //two electron indices
	count=0;
	int sw_count =0;

	for(auto s2=0; s2!=obs.size(); ++s2)
		for(auto s1=s2; s1!=obs.size(); ++s1)
			for(auto s4=0; s4!=obs.size(); ++s4)
				for(auto s3=s4; s3!=obs.size(); ++s3)  {
					if (Schwarz[s1*obs.size()+s2] * Schwarz[s3*obs.size()+s4]<1E-12 ) {
						sw_count++;
						continue;
					}
					eri_engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
					auto ints_shellset_eri = buf_vec_eri[0];

					auto bf1 = shell2bf[s1];     // first basis function in first shell
					auto n1 = obs[s1].size();     // number of basis functions in first shell
					auto bf2 = shell2bf[s2];     // first basis function in second shell
					auto n2 = obs[s2].size();     // number of basis functions in second shell
					auto bf3 = shell2bf[s3];     // first basis function in first shell
					auto n3 = obs[s3].size();     // number of basis functions in first shell
					auto bf4 = shell2bf[s4];     // first basis function in second shell
					auto n4 = obs[s4].size();     // number of basis functions in second shell
					if (ints_shellset_eri!=NULL) {



						for(auto f1=0; f1!=n1; ++f1)
							for(auto f2=0; f2!=n2; ++f2)
								for(auto f3=0; f3!=n3; ++f3)
									for(auto f4=0; f4!=n4; ++f4) {
										//Tmat_libint[(bf1+f1)*obs.nbf() + (bf2+f2)] = ints_shellset_t[f1*n2+f2]
										integral = ints_shellset_eri[f1*n2*n3*n4+f2*n3*n4+f3*n4+f4];
										if (((bf1+f1) >= (bf2+f2)) && ((bf3+f3) >= (bf4+f4)) && (((bf1+f1)*(bf1+f1+1)/2 +(bf2+f2))>=((bf3+f3)*(bf3+f3+1)/2+(bf4+f4))))
										{
                                            if (std::fabs(integral)>1.0E-12) {
												(*intval)[count] = integral;
												(*intnums)[count*4+0] = bf1+f1;
												(*intnums)[count*4+1] = bf2+f2;
												(*intnums)[count*4+2] = bf3+f3;
												(*intnums)[count*4+3] = bf4+f4;
												count++;
											}
										}
									}
					}
				}
				std::cout << "Schwarz count:" <<sw_count<< '\n';
    nrofint        = count;
	resort_integrals(*intnums,*intval,count,sortcount);

}


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
	std::cout << sorted << " integrals after search for perm_1\n";
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
	std::cout << sorted << " integrals after search for perm_15\n";
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
	std::cout << sorted << " integrals after search for perm_1234\n";
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
	std::cout << sorted << " integrals after search for perm_1256\n";
	sortcount[3] = sorted;
	std::cout << nrofint << " integrals in total\n";
	std::cout << "-------------------------------------------------------------------------------\n";

}
