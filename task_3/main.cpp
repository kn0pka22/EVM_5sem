#include <iostream>
#include "matr.hpp"


int main(int argc, char** argv) {

    int n=1;
	int m=0;
	int k=-1;
    int p;       //total_threads

	
	char* fname;


	
	if (argc<5 || argc>6){ printf("Please enter argc=5 or argc=6!\n"); return -1;} 
    if (((sscanf(argv[1], "%d", &n) != 1) || (sscanf(argv[2], "%d", &p) != 1) 
    || (sscanf(argv[3], "%d", &m) != 1) || (sscanf(argv[4], "%d", &k) != 1)|| (argc<5) || (argc>6) || (k<0) || (m<0) || (n<1) || (k>4))){ 
        
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
	if ((k==0) && (argc==5)){
		std::cout<<"Please enter the name of file!"<<std::endl;
        
		return -1; 
	}

	if (argc == 6 ) fname=argv[5];
    else fname=nullptr;
	//fname=argv[4];
	
	double* array_orig=(double*)malloc(n*n*sizeof(double));  
	
    //printf("%p - %ld\n",(void*)array_orig,n*n*sizeof(double));
    if (!array_orig) {
		std::cout<<"\nInfo:   not enough memory!"<<std::endl;
       
		return -1;
	}
	
	double* array=(double*)malloc(n*n*sizeof(double));  
    //printf("Not enough memory!222222222222222\n");
	if (!array) {
        free(array_orig);
		std::cout<<"\nInfo:   not enough memory!"<<std::endl;
        
		return -1;
	}
	double* array_inv=(double*)malloc(n*n*sizeof(double));
    //printf("Not enough memory!33333333333\n");
    if (!array_inv) {
        free(array_orig);
        free(array);
		std::cout<<"\nInfo:   not enough memory!"<<std::endl;
       
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
    if (!x_k) {
        free(array_orig);
        free(array);
        free(array_inv);
		std::cout<<"\nInfo:   not enough memory!"<<std::endl;
       
		return -1;
	}


//----------------------------for 3 task--------------------


    ARGS* args; //array of arguments for created tasks
    pthread_t* threads; //array of created task identifiers
    double* ArrayForNorm;

    if (!(args=(ARGS*)malloc(p*sizeof(ARGS)))){
        std::cout<<"\nInfo:   not enough memory"<<std::endl;
        free(array_orig);
        free(array);
        free(array_inv);
        free(x_k);
        return -2;
    }

    if (!(threads=(pthread_t*)malloc(p*sizeof(pthread_t)))){
        std::cout<<"\nInfo:   not enough memory"<<std::endl;
        free(array_orig);
        free(array);
        free(array_inv);
        free(x_k);
        free(args);
        return -2;
    }
    if (!(ArrayForNorm=(double*)malloc(p*sizeof(double)))){
        std::cout<<"\nInfo:   not enough memory"<<std::endl;
        free(array_orig);
        free(array);
        free(array_inv);
        free(x_k);
        free(args);
        free(threads);
        return -2;
    }



    for (int i = 0; i < p; i++) {
        args[i].matrix = array;
        args[i].inv = array_inv;
        args[i].x_k = x_k;
        args[i].n = n;
        args[i].thread_num = i;
        args[i].total_threads = p;
        args[i].count=p;
        args[i].flag_error=1;
        args[i].matrix_orig=array_orig;
    }

    //double time = get_full_time();


    //p=total threads

    
    for (int i = 0; i < p; ++i) {
        if (pthread_create(threads + i, 0, effective_method, args + i)) {
            std::cout <<"\nInfo:   the thread " << i << " has not been created" <<std::endl;
            free(threads);
            free(args);
            free(x_k);
            free(array);
            free(array_orig);
            free(array_inv);
            free(ArrayForNorm);
            return - 3;
        }
    }

    for (int i = 0; i < p; ++i) {
        if (pthread_join(threads[i], 0)) {
            std::cout << "\nInfo:   the thread" << i << " has not been joined!" << std::endl;
            free(threads);
            free(args);
            free(x_k);
            free(array);
            free(array_orig);
            free(array_inv);
            free(ArrayForNorm);
            return - 3;
        }
    }
    
	//time = get_full_time() - time;
	
	double time=0.0;
    bool err=true;
	for (int i=0;i<p;++i){
        if (args[i].flag_error==0) {std::cout<<"\nInfo:   an error occurred while working with the method\n"; err=0;}
        if (args[i].time>time) time=args[i].time;
    }
    

    for (int i = 0; i < p; ++i) {
        if (pthread_create(threads + i, 0, residual, args + i)) {
            std::cout <<"\nInfo:   the thread " << i << " has not been created" <<std::endl;
            free(threads);
            free(args);
            free(x_k);
            free(array);
            free(array_orig);
            free(array_inv);
            free(ArrayForNorm);
            return - 3;
        }
    }

    for (int i = 0; i < p; ++i) {
        if (pthread_join(threads[i], 0)) {
            std::cout << "\nInfo:   the thread" << i << " has not been joined!" << std::endl;
            free(threads);
            free(args);
            free(x_k);
            free(array);
            free(array_orig);
            free(array_inv);
            free(ArrayForNorm);
            return - 3;
        }
    }
    
     //gauss_back_run(array,array_inv, n);
    //std::cout<<"--------------------------Your  matrix   a   ------------------------\n"<<std::endl;
    //print_matrix_spv(array,n,n, m);
    double nor=0;
	if (err){
        std::cout<<"--------------------------Your inverse matrix------------------------\n"<<std::endl;
        print_matrix_spv_row(array_inv,n,n, m);
        
        double time2 = get_full_time();         
       // nor=norma(n,array_orig,array_inv);
       nor=args->nor;
        time2 = get_full_time() - time2;
        std::cout<<"\nInfo:   duration for norm = "<<time2<<std::endl;
    }
        
    printf("\n %s : residual = %e elapsed = %.2f s = %d n = %d m = %d p = %d\n", argv[0], nor, time, k, n, m, p);
    
    free(threads);
    free(args);
    free(x_k);
	free(array);
	free(array_orig);
	free(array_inv);
    free(ArrayForNorm);
	return 0;

}
