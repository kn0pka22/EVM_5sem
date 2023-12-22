#include <iostream>
#include <cmath>
#include <ctime>  

typedef struct{
    double* matrix;
    double* x_k; 
    double* inv; 
    int n; // matrix size 
    int thread_num; // thread number
    int total_threads; 
    double time;
    int count = 0;  //to count how many times the task was called
    int flag_error=1;  // to detect errors
    double* ArrayForNorm;
} ARGS;

void print_matrix(double* mas,int n);

double matrix_norm(int n, double *a);

int is_zero_matrix(int n, double *a);

void print_matrix_spv(double* mas,int n1,int n2,int m);

double norma(int n, double* ar, double* inv);

double norm_vec(int n,double* vec);

double formula(int k,int n, int ii,int jj);

int input_matr(int n,double* arr, int k,const char* fname);

//void ineffective_method(double* a,double* inv, int n);

int gauss_back_run(double* a,double* inv, int n);

int effective_method(double* a,double* inv, int n,double* x_k,int thread_num, int total_threads,int& flag);

void synchronize(int total_threads);

void* effective_method(void* pa);

//int effective_method(double* a,double* inv, int n,double* x_k);

double get_full_time();

void print_matrix_spv_row(double* mas,int n,int m,int p);

void* residual(void* pa);

