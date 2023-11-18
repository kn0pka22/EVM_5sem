#include <iostream>
#include <cmath>
#include <ctime>  

void print_matrix(double* mas,int n);

double matrix_norm(int n, double *a);

int is_zero_matrix(int n, double *a);

void print_matrix_spv(double* mas,int n1,int n2,int m);

double norm(int n, double* ar, double* inv);

double norm_vec(int n,double* vec);

double formula(int k,int n, int ii,int jj);

int input_matr(int n,double* arr, int k,const char* fname);

void ineffective_method(double* a,double* inv, int n);

int gauss_back_run(double* a,double* inv, int n);

int effective_method(double* a,double* inv, int n,double* x_k);
