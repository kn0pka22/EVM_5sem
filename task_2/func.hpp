#include <iostream>
#include <cmath>

void print_matrix_spv(double* mas,int n,int m,int p);

double formula(int k,int n, int ii,int jj);

double matrix_norm(int n, double *a);

int input_matr(int n,double* arr, int k,const char* fname);

int to_tridiag_form(double *a, double* y,double* x_k, double* z,int n);

int to_tridiag_form(double *a, double* y,double* x_k, double* z,double* matrB,int n);

int SymMatr(double* a, int n);  

int recursion(double* matr, int n, double eps,double left, double right, double* lambda,int index);

int n_(double * matr, double lmbd, int n);

double SearchEps();