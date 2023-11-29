#include <iostream>
#include <cmath>

void print_matrix_spv(double* mas,int n,int m,int p);

double formula(int k,int n, int ii,int jj);

double matrix_norm(int n, double *a);

int input_matr(int n,double* arr, int k,const char* fname);

int to_tridiag_form(double *a, double* y,double* x_k, double* z,int n);

int sym_matr(double* a, int n);

int sign_changes(double* a,int n,double lmbd);

double search_m_eps();

int search_values(int n,double* a,double* lmbd_values,double eps);


double discrepancy_inv1(int n,double* a,double* lmbd_values);

double discrepancy_inv2(int n,double* a,double* lmbd_values);


