#include <iostream>
#include "matr.hpp"


int main(int argc, char** argv) {

    int n=1;
	int m=0;
	int k=-1;
	int w;

	//double* array;
	//double* array_inv;
	//double* array_orig;
	
	 char* fname;
	
	if (argc<4 || argc>5){ printf("Please enter argc=4 or argc=5!\n"); return -1;} 
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 
    || sscanf(argv[3], "%d", &k) != 1 || (argc<4) || (argc>5) || (k<0) || (m<0) || (n<1) || (k>4)){ 
        
        printf("Invalid input!\n \
        * n – matrix dimension, \n \
        * m – the number of output results in the matrix, \n \
        * k – specifies the number of the formula for preparing matrices (-1< k <5), should be equal to 0 when input matrix from file \n \
        * filename – the name of the file from which the matrix should be read. This argument is missing if k! = 0.\n\n\
        Please enter: ./a.out n m k(>0)   \n\
        or    enter: ./a.out n m k(=0) filename  \n");

        return -1;
    }
    
	if (n<m) m=n;
	if (k==0 && argc==4){
		std::cout<<"Please enter the name of file!"<<std::endl;
        
		return -1; 
	}

	if (argc == 5 ) fname=argv[4];
    else fname=nullptr;
	//fname=argv[4];
	
	double* array_orig=(double*)malloc(n*n*sizeof(double));  
	
    //printf("%p - %ld\n",(void*)array_orig,n*n*sizeof(double));
    if (!array_orig) {
		printf("Not enough memory!\n");
       
		return -1;
	}
	
	double* array=(double*)malloc(n*n*sizeof(double));  
    //printf("Not enough memory!222222222222222\n");
	if (!array) {
        free(array_orig);
		printf("Not enough memory!\n");
        
		return -1;
	}
	double* array_inv=(double*)malloc(n*n*sizeof(double));
    //printf("Not enough memory!33333333333\n");
    if (!array_inv) {
        free(array_orig);
        free(array);
		printf("Not enough memory!\n");
       
		return -1;
	}

    //std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	for (int i=0;i<n;++i){
		for (int j=0;j<n;++j){
			if (i==j) array_inv[i*n+j]=1.0;
			else array_inv[i*n+j]=0.0;
		}
	}
	
	
	
		
	
	if (input_matr(n,array_orig, k,fname)!=1){
        //printf("\n");
        free(array_orig);
        free(array);
        free(array_inv);
       
        return -1;
    }
    else{std::cout<<"\nInfo:   Your input is correct\n"<<std::endl;}
    
	//for (int i=0;i<n*n;++i) std::cout<<array_orig[i]<<std::endl;

	
	for (int i=0;i<n;++i){
		for (int j=0;j<n;++j){
			array[i*n+j]=array_orig[i*n+j];
			} 
	}
	//if (w<0) free(array_inv);
 

	std::cout<<"-------------------------Your original matrix------------------------\n"<<std::endl;
	//print_matrix(array,m);
	print_matrix_spv(array,n,n, m);

	std::cout<<"-------------In the process of calculating the inverse matrix--------\n"<<std::endl;
	
	//ineffective_method(array,array_inv,n);
   
    double* x_k=(double*)malloc((n)*sizeof(double)); 

	clock_t start=clock();
	w=effective_method(array,array_inv,n,x_k);
	clock_t end=clock();

    double duration =(double)(end-start)/CLOCKS_PER_SEC;
    
	//print_matrix(array,m);
	
	//print_matrix(array_inv,m);
    if (w==0){

		std::cout<<"--------------------------Your inverse matrix------------------------\n"<<std::endl;
        print_matrix_spv(array_inv,n,n, m);
        std::cout<<"~~~~~~~Duration~~~~~~: \n"<<duration<<std::endl;
	    std::cout<<"~~~~~~~~~NORMA~~~~~~~: \n"<<norm(n,array_orig,array_inv)<<std::endl;
        
    
    }
    //else std::cout<<"Info:   Matrix singular! (0..0)"<<std::endl;
    
    free(x_k);
	free(array);
	free(array_orig);
	free(array_inv);
	return 0;

}
