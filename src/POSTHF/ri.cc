#include "ri.h"



void calculate_ri(libint2::BasisSet &obs,libint2::BasisSet &dfbs, double* BPQ)
{

    auto shell2bf_ob = obs.shell2bf();
    auto shell2bf_df = dfbs.shell2bf();

    int naux_2 = dfbs.nbf();
    int nroao  = obs.nbf();

    double* J        = new double[naux_2*naux_2];
    double* eigval   = new double[naux_2];
    double* eigvec   = new double[naux_2*naux_2];
    double* eigvec_c = new double[naux_2*naux_2];
    double* B        = new double[naux_2*nroao*nroao];

    libint2::Engine eri2_engine(libint2::Operator::coulomb, dfbs.max_nprim(), dfbs.max_l());
    eri2_engine.set_braket(libint2::BraKet::xs_xs);

    const auto& buf_vec_eri_J = eri2_engine.results();
    memset(J,0,naux_2*naux_2*sizeof(double));

    for(auto s1=0; s1!=dfbs.size(); ++s1)
        for(auto s2=0; s2!=dfbs.size(); ++s2) {
            eri2_engine.compute(dfbs[s1],dfbs[s2]);
            auto ints_shellset_eri = buf_vec_eri_J[0];
            if (ints_shellset_eri!=NULL) {

                auto bf1 = shell2bf_df[s1]; // first basis function in first shell
                auto n1 = dfbs[s1].size(); // number of basis functions in first shell

                auto bf2 = shell2bf_df[s2]; // first basis function in second shell
                auto n2 = dfbs[s2].size(); // number of basis functions in second shell

                for(auto f1=0; f1!=n1; ++f1)
                    for(auto f2=0; f2!=n2; ++f2)
                        J[(bf1+f1)*naux_2 + (bf2+f2)] = ints_shellset_eri[f1*n2+f2];
            }
        }


    std::memcpy(eigvec,J,naux_2*naux_2*sizeof(double));
    LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',naux_2,eigvec,naux_2,eigval);
    std::memcpy(eigvec_c,eigvec,naux_2*naux_2*sizeof(double));



    double max_eig = eigval[naux_2-1];
    double tol = 1.0E-10;
    for (int i = 0; i<naux_2; i++) {
        if(eigval[i]/max_eig < tol || eigval[i] <0)
            eigval[i] =0.0;
        else{
            eigval[i] = 1 / std::sqrt(eigval[i]);
        }

        cblas_dscal(naux_2,eigval[i],&eigvec[i],naux_2);
    }

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,naux_2,naux_2,
                naux_2,1.0,eigvec_c,naux_2,eigvec,naux_2,0.0,J,naux_2);
    libint2::Engine eri3_engine(libint2::Operator::coulomb,
                                std::max(dfbs.max_nprim(),obs.max_nprim()),
                                std::max(dfbs.max_l(),obs.max_l()));
    memset(B,0,naux_2*nroao*nroao*sizeof(double));
    eri3_engine.set_braket(libint2::BraKet::xs_xx);
    const auto& buf_vec_eri3 = eri3_engine.results();

    for(auto s1=0; s1!=dfbs.size(); ++s1) {
        for (auto s2=0; s2!= obs.size(); ++s2) {
            for (auto s3=0; s3!=obs.size(); ++s3) {
                eri3_engine.compute(dfbs[s1],obs[s2],obs[s3]);
                auto ints_shellset_eri3 = buf_vec_eri3[0];

                if (ints_shellset_eri3!=NULL) {
                    auto bf1 = shell2bf_df[s1]; // first basis function in first shell
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

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,naux_2,nroao*nroao,naux_2,1.0,J,naux_2,B,nroao*nroao,0.0,
                BPQ,nroao*nroao);

    delete[] eigval;
    delete[] eigvec_c;
    delete[] eigvec;
    delete[] J;
    delete[] B;
}


void read_transform_ri(std::string prefix,  //prefix to find the file
                       int nroe,            //nr of electrons
                       int nroao,           //nr of bsf in aobasis
                       int naux_2,          //nr of aux basis functions
                       double* MOs,         //orbitals for transformation
                       double* Bia)        //output - transformed array containing mointegrals
{

    int nocc = nroe/2;
    int nvir = nroao-nroe/2;
    double* BPQ = new double[naux_2*nroao*nroao];
    double* BiQ = new double[naux_2*nocc*nroao];

    memset(BiQ,0,sizeof(double)*naux_2*nocc*nroao);
    memset(Bia,0,sizeof(double)*naux_2*nocc*nvir);

    std::ifstream datf;
    datf.open(prefix+".rimp2");
    datf.read((char*) BPQ, naux_2*nroao*nroao*sizeof(double));
    datf.close();

    for (int Q=0; Q<naux_2; Q++)
        for (int i=0; i<nocc; i++)
            for (int mu=0; mu<nroao; mu++)
                for (int nu=0; nu<nroao; nu++)
                    BiQ[Q*nocc*nroao + i*nroao + nu] +=  MOs[i*nroao+mu]*BPQ[Q*nroao*nroao + mu*nroao + nu];
    for (int Q=0; Q<naux_2; Q++)
        for (int a=0; a<nvir; a++)
            for (int i =0; i<nocc; i++)
                for (int nu=0; nu<nroao; nu++)
                    Bia[Q*nocc*nvir+i*nvir+a]+= MOs[(a+nroe/2)*nroao + nu]*BiQ[Q*nocc*nroao + i*nroao + nu ];


    delete[] BiQ;
    delete[] BPQ;

}

void transform_ri(int nroe,            //nr of electrons
                  int nroao,           //nr of bsf in aobasis
                  int naux_2,          //nr of aux basis functions
                  double* MOs,         //orbitals for transformation
                  double* BPQ,          //input - calculated b^Q_pq
                  double* Bia)         //output - transformed array containing mointegrals)
{
    int nocc = nroe/2;
    int nvir = nroao-nroe/2;
    double* BiQ = new double[naux_2*nocc*nroao];
    memset(BiQ,0,sizeof(double)*naux_2*nocc*nroao);
    memset(Bia,0,sizeof(double)*naux_2*nocc*nvir);
    for (int Q=0; Q<naux_2; Q++)
        for (int i=0; i<nocc; i++)
            for (int mu=0; mu<nroao; mu++)
                for (int nu=0; nu<nroao; nu++)
                    BiQ[Q*nocc*nroao + i*nroao + nu] +=  MOs[i*nroao+mu]*BPQ[Q*nroao*nroao + mu*nroao + nu];
    for (int Q=0; Q<naux_2; Q++)
        for (int a=0; a<nvir; a++)
            for (int i =0; i<nocc; i++)
                for (int nu=0; nu<nroao; nu++)
                    Bia[Q*nocc*nvir+i*nvir+a]+= MOs[(a+nroe/2)*nroao + nu]*BiQ[Q*nocc*nroao + i*nroao + nu ];
    delete[] BiQ;
}
