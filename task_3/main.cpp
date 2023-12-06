#include <iostream>
#include "matr.hpp"





int main(int argc, char** argv) {

    int n=1;
	
	int m=0; //r
	int k=-1;
	int w;
	int thread_count=1;  // p
	long double time;

	//double* array;
	//double* array_inv;
	//double* array_orig;
	
	char* fname;
	
	if (argc<5 || argc>6){ printf("Please enter argc=5 or argc=6!\n"); return -1;} 
    if (((sscanf(argv[1], "%d", &n) != 1) || (sscanf(argv[2], "%d", &thread_count) != 1 )
    || (sscanf(argv[3], "%d", &m) != 1) ||  (sscanf(argv[4], "%d", &k) != 1) || (argc<5) || (argc>6) || (k<0) || (m<0) || (n<1) ||   (thread_count<1) || (k>4))){ 
        
        printf("Invalid input!\n \
        * n – matrix dimension, \n \
		* p - number of threads used, \n \
        * m – the number of output results in the matrix, \n \
        * k – specifies the number of the formula for preparing matrices (-1< k <5), should be equal to 0 when input matrix from file \n \
        * filename – the name of the file from which the matrix should be read. This argument is missing if k! = 0.\n\n\
        Please enter: ./a.out n m k(>0)   \n\
        or    enter: ./a.out n m k(=0) filename  \n");

        return -1;
    }
    
	if (n<m) m=n;
	if (k==0 && argc==5){
		std::cout<<"Please enter the name of file!"<<std::endl;
        
		return -1; 
	}

	if (argc == 6 ) fname=argv[5];
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

//-----------------------------
	
    

	pthread_t* threads=(pthread_t*)malloc((thread_count)*sizeof(pthread_t)); 
	

	if (!threads){
		free(array_orig);
        free(array);
        free(array_inv);
		free(x_k);
		return -2;
	}

	ARGS* args=(ARGS*)malloc((thread_count)*sizeof(ARGS));

	if (!args){
		free(array_orig);
        free(array);
        free(array_inv);
		free(x_k);
		free(threads);
		return -2;
	}

	for (int i = 0; i < thread_count; ++i) {
        args[i].a = array;
		args[i].x_k = x_k;
		args[i].n = n;
		args[i].thread_num = i;
		args[i].total_threads = thread_count;
		args[i].k = k;
    }

//int pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*start_routine)(void *), void *arg)

/*
The thread parameter is of type pthread_t. The created thread will always be accessible with this reference.
The attr parameter lets you specify custom behavior. You can use a series of thread-specific functions starting with pthread_attr_ to set up this value. Possible customizations are the scheduling policy, stack size, and detach policy.
start_routine specifies the function that the thread will run.
arg represents a generic data structure passed to the function by the thread.
*/


//std::cout<<"thread_count = "<<thread_count<<std::endl;
	for (int i = 0; i < thread_count; i++){
		if (pthread_create(threads + i, 0, effective_method, args + i)){
			std::cout << "the thread with number " << i << " has not been created" << std::endl;
			free(threads);
			free(args);
			free(x_k);
			free(array);
			free(array_orig);
			free(array_inv);
		}
	}

	for (int i = 0; i < thread_count; i++) {
		
			if (pthread_join(threads[i], 0)) {
				std::cout << "the thread with number " << i << " has not been started" << std::endl;
				free(threads);
				free(args);
				free(x_k);
				free(array);
				free(array_orig);
				free(array_inv);
			}
			
	}

//----------------------------

	//clock_t start=clock();
	//w=effective_method(array,array_inv,n,x_k);
	//clock_t end=clock();

    double duration =1;//=(double)(end-start)/CLOCKS_PER_SEC;
    
	//print_matrix(array,m);
	
	//print_matrix(array_inv,m);
	printf("%s : residual = %e elapsed = %.2f s = %d n = %d m = %d p = %d\n", argv[0], norm(n,array_orig,array_inv), duration, k, n, m, thread_count);

    /*if (w==0){

		std::cout<<"--------------------------Your inverse matrix------------------------\n"<<std::endl;
        print_matrix_spv(array_inv,n,n, m);
        std::cout<<"Info:   duration                             : "<<duration<<std::endl;
        std::cout<<"Info:   norma                                : "<<norm(n,array_orig,array_inv)<<std::endl;
    
    
    }*/
    //else std::cout<<"Info:   Matrix singular! (0..0)"<<std::endl;
    

	free(threads);
	free(args);
    free(x_k);
	free(array);
	free(array_orig);
	free(array_inv);
	return 0;

}
