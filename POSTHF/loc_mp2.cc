#include <iostream>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <cmath>

void run_non_canonical_mp2(int nroao, int nroe, double* prec_ints, double* FMo)
{
    std::cout << "\nCalculating semi-canonical amplitudes ... ";

    int nocc = nroe/2;
    int nvir = nroao - nocc;
    double* T_ijab   = new double[nocc*nocc*nvir*nvir];

    int Ti,Tj,Ta,Tb = 0;
    int Ta_step = nvir;
    int Tj_step = Ta_step * nvir;
    int Ti_step = Tj_step * nocc;

    long long int kstep = nroao;
    long long int jstep = nroao*kstep;
    long long int istep = nroao*jstep;

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
}
