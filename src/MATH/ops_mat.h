#include <complex>

void diag_mat(int nroao, double* mat, double* vals, double* vecs);
void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);         //for large matrices
void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir); //for large matrices
void pmv(double* mat, double* vi, double* vo, int nroao);
