
#include "func.hpp"



int main(int argc, char** argv){
	int n=1;
	int m=0;
	int k=-1;
	double eps=100;
	char* fname;
	
	if (argc<5 || argc>6){ std::cout<<"Please enter argc=5 or argc=6!\n"; return -1;} 
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%lf", &eps) != 1
    || sscanf(argv[4], "%d", &k) != 1 || (argc<5) || (argc>6) || (k<0) || (m<0) || (n<1) || (k>5) || (eps<0.0)){ 
        
        std::cout<<"Invalid input!\n \
        * n вЂ“ matrix dimension, \n \
        * m вЂ“ the number of output results in the matrix, \n \
        * k вЂ“ specifies the number of the formula for preparing matrices (-1< k <5), should be equal to 0 when input matrix from file \n \
        * filename вЂ“ the name of the file from which the matrix should be read. This argument is missing if k! = 0.\n\n\
        Please enter: ./a.out n m k(>0)   \n\
        or    enter: ./a.out n m k(=0) filename  \n";

        return -1;
    }
    
	if (n<m) m=n;
	if (k==0 && argc==5){
		std::cout<<"Please enter the name of file!"<<std::endl;
    	return -1; 
	}

	if (argc == 6 ) fname=argv[5];
    else fname=nullptr;
	
	double* array_orig=(double*)malloc(n*n*sizeof(double));  
	
    //printf("%p - %ld\n",(void*)array_orig,n*n*sizeof(double));
    if (!array_orig) {
		printf("Not enough memory!\n");
		return -1;
	}
	
	double* array=(double*)malloc(n*n*sizeof(double)); 
	if (!array) {
        free(array_orig);
		printf("Not enough memory!\n");
		return -1;
	}
		
	
	if (input_matr(n,array_orig, k,fname)!=1){
        free(array_orig);
        free(array); 
        return -1;
    }
    else{std::cout<<"\nInfo:   Your input is correct\n"<<std::endl;}
    

	for (int i=0;i<n;++i){
		for (int j=0;j<n;++j){
			array[i*n+j]=array_orig[i*n+j];
			} 
	}
	

	std::cout<<"----------------------------------Your original matrix--------------------------------\n"<<std::endl;
	print_matrix_spv(array,n,n, m);

	
	double* x_k=(double*)malloc(n*sizeof(double)); 
	if (!x_k){
		free(array_orig);
		free(array);
		std::cout<<"Info:   Not enough memory!\n";
		return -1;
	}

	double* y=(double*)malloc(n*sizeof(double)); 
	if (!y){
		free(array_orig);
		free(array);
		free(x_k);
		std::cout<<"Info:   Not enough memory!\n";
		return -1;
	}

	
	double* z=(double*)malloc(n*sizeof(double)); 
	if (!z){
		free(array_orig);
		free(array);
		free(y);
		free(x_k);
		std::cout<<"Info:   Not enough memory!\n";
		return -1;
	}


	double* lmbd_values=(double*)malloc(n*n*sizeof(double)); 
	if (!lmbd_values){
		free(array_orig);
		free(array);
		free(y);
		free(x_k);
		free(z);
		std::cout<<"Info:   Not enough memory!\n";
		return -1;
	}

	
	if (sym_matr(array,n)<0){
		free(array_orig);
		free(array);
		free(y);
		free(x_k);
		free(z);
		std::cout<<"Info:   the matrix is not symmetrical!\n";
		return -1;
	}



	double t1 = clock();
	to_tridiag_form(array, y,x_k,z,n);
	t1 = (clock() - t1) / CLOCKS_PER_SEC;


	
	double t2 = clock();
	search_values(n,array,lmbd_values, eps);
	t2 = (clock() - t2) / CLOCKS_PER_SEC;

	std::cout<<"--------------------------------------- tridiag form ---------------------------------\n"<<std::endl;
	print_matrix_spv(array,n,n,m);
	//int num=sign_changes(array,n,0);
	//std::cout<<"Info:   number of negative eigenvalues = "<<num<<std::endl;
	//
	std::cout<<"Info:   duration for tridiag method          : "<<t1<<std::endl;
	std::cout<<"Info:   duration of searching for eigenvalues: "<<t2<<std::endl;


	//double duration =(double)(t1+t2)/CLOCKS_PER_SEC;
	std::cout<<"Info:   duration                             : "<<(t1+t2)<<std::endl;


	double tmp=lmbd_values[0];
	std::cout<<"---------------------------------- Matrix eigenvalues ------------------------------- \n"<<std::endl;
	for (int i=0;i<n;++i){
		if ((lmbd_values[i]-tmp)<eps) std::cout<<i+1<<")"<<lmbd_values[i]<<"   ";
		else std::cout<<std::endl<<i+1<<")"<<lmbd_values[i]<<"   ";
		tmp=lmbd_values[i];
	}
	std::cout<<"\n\n"<<std::endl;
    
    //std::cout<<"HEREE!!!"<<std::endl;
	
  
        
	std::cout<<"Info:   discrepancy_inv1 = "<<discrepancy_inv1(n,array,lmbd_values)<<std::endl;
	std::cout<<"Info:   discrepancy_inv2 = "<<discrepancy_inv2(n,array,lmbd_values)<<std::endl;
	//---------------------------not for task------------
	
	/*
	char* A=new char[n];
	char* B=new char[n];
	char* C=new char[n];
	
	for (int i=0;i<n;++i ) B[i] = rand();
	for (int i=0;i<n;++i ) C[i] =rand(); 

	std::cout<<"\n array B:  "; for (int i=0;i<n;++i ) std::cout<<B[i]<<"  ";
	std::cout<<"\n array C:  "; for (int i=0;i<n;++i ) std::cout<<C[i]<<"  ";

	double t_foo1 = clock();
	foo(A, B);
	t_foo1 = (clock() - t_foo1) / CLOCKS_PER_SEC;

	double t_foo2 = clock();
	foo2(A, C);
	t_foo2 = (clock() - t_foo2) / CLOCKS_PER_SEC;
	std::cout<<"\n array A:  "; for (int i=0;i<n;++i ) std::cout<<A[i]<<"  ";
	std::cout<<std::endl;
	std::cout<<"Duration for foo1: "<<t_foo1<<" and for foo2: "<<t_foo2<<std::endl;

	
	delete[] A;
	delete[] B;
	delete[] C;
	
*/

	free(array);
	free(array_orig);
	free(y);
	free(z);
	free(x_k);
	free(lmbd_values);
	


	return 0;
}
