#include "func.hpp"



int main(int argc, char** argv){
	int n=1;
	int m=0;
	int k=-1;
	double eps=100;
	char* fname;
	
	if (argc<5 || argc>6){ std::cout<<"Please enter argc=5 or argc=6!\n"; return -1;} 
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%lf", &eps) != 1
    || sscanf(argv[4], "%d", &k) != 1 || (argc<5) || (argc>6) || (k<0) || (m<0) || (n<1) || (k>4) || (eps<0.0)){ 
        
        std::cout<<"Invalid input!\n \
        * n – matrix dimension, \n \
        * m – the number of output results in the matrix, \n \
        * k – specifies the number of the formula for preparing matrices (-1< k <5), should be equal to 0 when input matrix from file \n \
        * filename – the name of the file from which the matrix should be read. This argument is missing if k! = 0.\n\n\
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
		
	
	if (input_matr(n,array_orig, k,fname)!=1 || SymMatr(array_orig, n)!=1){
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
	

	std::cout<<"-------------------------Your original matrix------------------------\n"<<std::endl;
	print_matrix_spv(array,n,n, m);

	
	double* x_k=(double*)malloc(n*sizeof(double)); 
	if (!x_k){
		free(array_orig);
		free(array);
		std::cout<<"Not enough memory!\n";
		return -1;
	}

	double* y=(double*)malloc(n*sizeof(double)); 
	if (!y){
		free(array_orig);
		free(array);
		free(x_k);
		std::cout<<"Not enough memory!\n";
		return -1;
	}

	
	double* z=(double*)malloc(n*sizeof(double)); 
	if (!z){
		free(array_orig);
		free(array);
		free(y);
		free(x_k);
		std::cout<<"Not enough memory!\n";
		return -1;
	}


	double* matrB=(double*)malloc(n*n*sizeof(double)); 
	if (!matrB){
		free(array_orig);
		free(array);
		free(y);
		free(x_k);
		free(z);
		std::cout<<"Not enough memory!\n";
		return -1;
	}

	double* lmbd_values=(double*)malloc(n*n*sizeof(double)); 
	if (!matrB){
		free(array_orig);
		free(array);
		free(y);
		free(x_k);
		free(z);
		free(matrB);
		std::cout<<"Not enough memory!\n";
		return -1;
	}

	for (int i=0;i<n;++i){
		for (int j=0;j<n;++j){
			matrB[i*n+j]=array_orig[j*n+i];
		}   
	}

	double norma = matrix_norm(n,array);
	std::cout<<"-------------In the process of calculating --------\n"<<std::endl;
	
	int w = to_tridiag_form(array, y,x_k,z, matrB,n);
	//w=recursion(matrB, n, eps,(-1.0)*norma, norma, lmbd_values,0);
	//n_(matrB, 1, n);

	print_matrix_spv(matrB,n,n, m);
	//n_(array, 1, n);

	/*
	std::cout<<"values: ";
	for (int i=0;i<n;++i){ 
		std::cout<<lmbd_values[i]<<" ";
	}
	std::cout<<std::endl;

	*/
    
	free(array);
	free(array_orig);
	free(y);
	free(z);
	free(x_k);
	free(matrB);
	free(lmbd_values);
	


	return 0;
}
