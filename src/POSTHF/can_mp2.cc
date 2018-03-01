#include <iostream>
#include <cstring>
#include <iomanip>
#include <cblas.h>
#include "../IO/ops_io.h"
#define PWIDTH_L 12
#define PWIDTH_R 16



void run_canonical_mp2(systeminfo* sysinfo,pHF *postHF)             //Fock Matrix in MO Basis
{
    int nroe = sysinfo->nroe;
    int nroao = sysinfo->nroao;
    double* prec_ints = postHF->prec_ints;
    double* FMo = postHF->FMo;
    double EMP2 = 0.0;
    double EMP2_SS = 0.0;
    double EMP2_OS = 0.0;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;


    for (int i = 0; i < nroe/2; i++) {
      for (int j = 0; j < nroe/2; j++) {
        for (int a = nroe/2; a < nroao; a++) {
          for (int b = nroe/2;  b < nroao; b++) {
            EMP2_OS += - (prec_ints[i*istep + a*jstep + j*kstep + b]*(prec_ints[i*istep + a*jstep + j*kstep + b]))/ (FMo[a*nroao +a]+FMo[b*nroao+b]-FMo[i*nroao+i]-FMo[j*nroao+j]);
            EMP2_SS += - (prec_ints[i*istep + a*jstep + j*kstep + b] - prec_ints[i*istep + b*jstep + j*kstep + a])*prec_ints[i*istep + a*jstep + j*kstep + b]/(FMo[a*nroao+a]+FMo[b*nroao+b]-FMo[i*nroao +i]-FMo[j*nroao + j]);
          }
        }
      }
    }
  EMP2 = EMP2_OS + EMP2_SS;
  std::cout<<"\nMP2-Results:\n";
  std::cout<<"------------\n";
  std::cout <<std::setw( PWIDTH_L ) << "EMP2_SS:" <<std::setw( PWIDTH_R )<<EMP2_SS<< '\n';
  std::cout <<std::setw( PWIDTH_L ) << "EMP2_OS:" <<std::setw( PWIDTH_R )<<EMP2_OS<< '\n';
  std::cout <<std::setw( PWIDTH_L ) << "EMP2:"    <<std::setw( PWIDTH_R )<<EMP2<< '\n';
  std::cout<<"\nSCS-MP2-Results:\n";
  std::cout<<"------------\n";
  std::cout <<std::setw( PWIDTH_L ) << "EMP2_SS:" <<std::setw( PWIDTH_R )<<EMP2_SS/3.0 << '\n';
  std::cout <<std::setw( PWIDTH_L ) << "EMP2_OS:" <<std::setw( PWIDTH_R )<<EMP2_OS*6.0/5.0 << '\n';
  std::cout <<std::setw( PWIDTH_L ) << "EMP2(SCS):"    <<std::setw( PWIDTH_R )<<EMP2_SS/3. +EMP2_OS*6.0/5.0 << '\n';


}



void run_canonical_mp2_ri(int nroe,         //Number of Electrons
                          int nroao,        //Nr of bsf in orbital basis
                          int naux_2,       //Nr of bsf in ri-basis
                          double* Bia,      //tranformed intergals
                          double* FMo)      //Fock Matrix in MO Basis
{

    int nocc = nroe/2;
    int nvir = nroao-nroe/2;

    double iajb = 0.0;
    double ibja = 0.0;
    double eps_i = 0.0;
    double eps_j = 0.0;
    double eps_a = 0.0;
    double eps_b = 0.0;
    double EMP2,EMP2_SS,EMP2_OS;

    EMP2=EMP2_SS=EMP2_OS=0.0;

    double* locbuff = new double[4*nvir*nvir];
    int inc=0;

    double* M       = &locbuff[inc]; inc+=nvir*nvir;
    double* MD      = &locbuff[inc]; inc+=nvir*nvir;
    double* X       = &locbuff[inc]; inc+=nvir*nvir;
    double* tmp     = &locbuff[inc]; inc+=nvir*nvir;

    for (int i = 0; i < nroe/2; i++) {
        eps_i = FMo[i*nroao+i];

        for (int j = 0; j < nroe/2; j++) {
              eps_j = FMo[j*nroao+j];

     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nvir,nvir,naux_2,1.0,&Bia[i*nvir],nvir*nocc,&Bia[j*nvir],nvir*nocc,0.0,M,nvir);
     for (int a = 0; a < nvir; a++)  {
         for (int b = 0;  b < nvir; b++){
             X[a*nvir+b] = 1/(FMo[(a+nocc)*nroao+(a+nocc)]+FMo[(b+nocc)*nroao+(b+nocc)]-eps_i-eps_j);
             MD[a*nvir+b] = M[a*nvir+b]-M[b*nvir+a];
        }
     }
     cblas_dsbmv(CblasRowMajor,CblasLower,nvir*nvir,0,1.0,M,1,X,1,0.0,tmp,1);
     EMP2_OS += - cblas_ddot(nvir*nvir,tmp,1,M,1);
     EMP2_SS += - cblas_ddot(nvir*nvir,tmp,1,MD,1);

     //cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,nvir,nvir,nvir,1.0,M,nvir,X,nvir,0.0,tmp2,nvir);
     //EMP2_OS += - cblas_ddot(nvir*nvir,tmp,1,M,1);
     //for (int a=0;a<nvir;a++)
     //     EMP2_OS += tmp2[i*nvir+i];
    /*
     for (int a = 0; a < nvir; a++)  {
        for (int b = 0;  b < nvir; b++){
            //EMP2_OS += - ((M[a*nvir+b] - M[b*nvir+a])*tmp[a*nvir+b]);
          }
        }
     */

     }
   }

    EMP2 = EMP2_OS + EMP2_SS;
    std::cout<<"\nRI-MP2-Results:\n";
    std::cout<<"---------------\n";
    std::cout <<std::setw( PWIDTH_L ) << "EMP2_SS:" <<std::setw( PWIDTH_R )<<EMP2_SS<< '\n';
    std::cout <<std::setw( PWIDTH_L ) << "EMP2_OS:" <<std::setw( PWIDTH_R )<<EMP2_OS<< '\n';
    std::cout <<std::setw( PWIDTH_L ) << "EMP2:"    <<std::setw( PWIDTH_R )<<EMP2<< '\n';

    std::cout<<"\nSCS-RI-MP2-Results:\n";
    std::cout<<"-------------------\n";
    std::cout <<std::setw( PWIDTH_L ) << "EMP2_SS:" <<std::setw( PWIDTH_R )<<EMP2_SS/3.0 << '\n';
    std::cout <<std::setw( PWIDTH_L ) << "EMP2_OS:" <<std::setw( PWIDTH_R )<<EMP2_OS*6.0/5.0 << '\n';
    std::cout <<std::setw( PWIDTH_L ) << "EMP2(SCS):"    <<std::setw( PWIDTH_R )<<EMP2_SS/3. +EMP2_OS*6.0/5.0 << '\n';


    delete[] locbuff;



}
