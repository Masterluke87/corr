#include <iostream>
#include <iomanip>
#define PWIDTH_L 12
#define PWIDTH_R 16

void run_canonical_mp2(int nroe,int nroao,double* prec_ints,double* FMo)
{
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
            EMP2_SS += - (prec_ints[i*istep + a*jstep + j*kstep + b]*(prec_ints[i*istep + a*jstep + j*kstep + b]))/ (FMo[a*nroao +a]+FMo[b*nroao+b]-FMo[i*nroao+i]-FMo[j*nroao+j]);
            EMP2_OS += - (prec_ints[i*istep + a*jstep + j*kstep + b] - prec_ints[i*istep + b*jstep + j*kstep + a])*prec_ints[i*istep + a*jstep + j*kstep + b]/(FMo[a*nroao+a]+FMo[b*nroao+b]-FMo[i*nroao +i]-FMo[j*nroao + j]);
          }
        }
      }
    }
  EMP2 = EMP2_OS + EMP2_SS;
  std::cout <<std::setw( PWIDTH_L ) << "EMP2_SS:" <<std::setw( PWIDTH_R )<<EMP2_SS<< '\n';
  std::cout <<std::setw( PWIDTH_L ) << "EMP2_OS:" <<std::setw( PWIDTH_R )<<EMP2_OS<< '\n';
  std::cout <<std::setw( PWIDTH_L ) << "EMP2:"    <<std::setw( PWIDTH_R )<<EMP2<< '\n';
}
