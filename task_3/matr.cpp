#include "matr.hpp"
#include <chrono>
#include <ctime>
#include <sys/time.h>
std::mutex mtx;

double get_full_time(){
	struct timeval buf;
	gettimeofday(&buf, 0);    //преобразуем время в секундах в сотые доли секунды
	return buf.tv_sec  + buf.tv_usec/ 1000000.0;
}




void print_matrix(double* mas,int n){

	if (mas==NULL) return;
	for(int i = 0;i<n; i++) {
		for(int j = 0;j<n;j++){
	    		printf(" %10.3e ",mas[j*n+i]);
		}
	    printf("\n");
    }
	printf("\n\n");
}

void print_matrix_spv(double* mas,int n,int m,int p){
    
    int i, j, n_max = (n > p ? p : n), m_max = (m > p ? p : m);
    for (i = 0; i < n_max; i ++) {
		printf("     ");
		for (j = 0; j < m_max; j ++)
			printf(" %10.3e", mas[j*m + i]);
		printf("\n");
	}
    printf("\n\n");
}

void print_matrix_spv_row(double* mas,int n,int m,int p){  //for a matrix stored row by row
    
    int i, j, n_max = (n > p ? p : n), m_max = (m > p ? p : m);
    for (i = 0; i < n_max; i ++) {
		printf("     ");
		for (j = 0; j < m_max; j ++)
			printf(" %10.3e", mas[i*m + j]);
		printf("\n");
	}
    printf("\n\n");
}

double formula(int k,int n, int ii,int jj){
	int p,q,i,j;
	i=ii+1;
	j=jj+1;
	if (i>j) {p=i;q=j;}
	else {p=j;q=i;}
	
	if     (k==1) return n-p+1;
	else if(k==2) return p; 
	else if(k==3) return p-q;
	else if(k==4) return 1.0/(p+q-1);
	return 0;
}




double norma(int n, double* ar, double* inv) {
	double ma = 0.0;
	double str_sum = 0.0;
	double sum=0.0;

	for (int i = 0; i < n; i++) {
		str_sum = 0.0;
		for (int j = 0; j < n; j++) {
			sum = 0.0;
			for (int k = 0; k < n; k++) {
				sum += (ar[k * n + j] * inv[k * n + i]);
			}
			if (i == j) { sum -= 1.0; }
			
			str_sum += fabs(sum);
		}
		if (fabs(str_sum) > ma) ma = str_sum;
	}
	return ma;
}


double norm_vec(int m,double* vec){
	double sum=0.0;
	for (int i=0;i<m;++i){
		sum+=vec[i]*vec[i];
	}
	return sqrt(sum);
}

double matrix_norm(int n, double *a){
	double sum=0;
    double max=-10000;
	double norm_m;

   	for(int i = 0; i < n; i++)
    { 
        for(int j = 0; j < n; j++)
            sum += fabs(a[j*n+i]); 
        if(max < sum)
            max = sum;
        sum = 0.0; 
    } 
    norm_m = max; 
    return norm_m;
}



int is_zero_matrix(int n, double *a){

	int cnt=0;
	int cnt2=-1;
	int cnt3=-1;
	int cnt4=0;
	double tmp,tmp3;
    double norm = matrix_norm(n,a);
    for (int i = 0; i < n-1; ++i) { 
		tmp=std::abs(a[i * n + i]);
		if  ((tmp< (1e-15*norm)) || ((a[i * n + i]-0.0)<1e-15*norm)) return -1;
    }
       // std::cout<<"HERE!"<<std::endl;
	tmp=a[(n-1)*n];
	tmp3=a[n-1];
	for (int i=0;i<n;++i){
		if ((std::abs(a[(n-1)*n+i])< 1e-15*norm) ) cnt++; 
		if  ((a[(n-1)*n+i]-tmp)<1e-15*norm) cnt2++; 
		if  ((a[i*(n-1)]-tmp3)<1e-15*norm) cnt3++; 
		if ((std::abs(a[i*n+n-1])< 1e-15*norm)) cnt4++; 
	}
	//std::cout<<"cnt="<<cnt<<" and cnt2="<<cnt2<<" and cnt3="<<cnt3<<" and cnt4="<<cnt4<<std::endl;
	if (cnt==n || cnt2==n || cnt3==n || cnt4==n) return -1;
	
    //std::cout<<"MATRIX SINGULAR!\n";
    return 1;
}



int input_matr(int n,double* arr, int k,const char* fname){
		double a;
		int q;
		char r;
        int cnt=0;
        
		if (k==0){
			FILE* f;
			f=fopen(fname,"r");
			if (!f){ 
				std::cout<<"\nInfo:   File not open(\n";
				//free(arr);
				return -2;
			}
			
			else{
				while(!feof(f)){
					q=fscanf(f,"%s",&r);
					//std::cout<<r<<std::endl;
					if((r >= 48 && r <= 57) || (r=='-'));
					else {
						std::cout<<"\nInfo:   unrecognized character in file!\n"; 
						return -1;
					}
				}
				fseek(f, 0, SEEK_SET);
				
				for (int i=0;i<n;++i){
					for (int j=0;j<n;++j){		
						q=fscanf(f,"%lf",&a);
						if (q!=1){
							
							std::cout<<"\nInfo:   not enough elements in the file..\n";
							fclose(f);
							//free(arr);
							return -2;
						}
						arr[j*n+i]=a;
						//std::cout<<"array["<<i*n+j<<"]="<<array[i*n+j]<<std::endl;
					}
				}
				fclose(f);
			}
			double norm = matrix_norm(n,arr);
			//if(is_zero_matrix(n, arr)==-1) return -1;
			for (int i=0;i<n*n;++i){ if (std::abs(arr[i]-0.0)<=1e-15*norm){ cnt++;}  }
            if (cnt==n*n) {std::cout<<"\nInfo:   Singular matrix!"<<std::endl; return -1;}
		}
		else{
			for (int i=0;i<n;++i){
				for (int j=0;j<n;++j){		
					arr[j*n+i]=formula(k,n,i,j);
				}
			}
		
		}
		
	return 1;
}



int gauss_back_run(double* a,double* inv, int n){
 //std::cout<<"Here!";
   
    double tmp;
    for (int k=n-1;k>=1;--k){
		for (int i=k-1;i>=0;--i){
			tmp=a[k*n+i];
            //std::cout<<"NORMA: "<<norm_m*1e-15<<std::endl;
            //if (fabs(a[k*n+k]))<((1e-15)*norm_m)){ std::cout<<"INCORRECT MATRIX!"<<std::endl; return -1;}
            //if (fabs(a[k*n + k])<=((1e-15)*norm_m)) {std::cout<<"INCORRECT MATRIX!"<<std::endl; return -1;}
            double MEGAVAR = -tmp/a[k*n+k];
			for (int j=n-1;j>-1;--j){
				a[j*n+i]+=MEGAVAR*a[j*n+k];
				inv[i*n+j]+=MEGAVAR * inv[k*n+j];
			}
			//std::cout<<"matr a: \n"; print_matrix_spv(a,n,n,5);
			//std::cout<<"matr inv: \n"; print_matrix_spv(inv,n,n,5);
		}
	}
	for (int i=0;i<n;++i){
		tmp=a[i*n+i];
        //std::cout<<(1e-15)*norm_m<<std::endl;
        //if (fabs(tmp)<=((1e-15)*norm_m)) {std::cout<<"SINGULAR MATRIX!"<<std::endl; return -1;}
		for (int j=0;j<n;++j){
			a[j*n+i]=a[j*n+i]/tmp;
			inv[i*n+j]=inv[i*n+j]/tmp;
		}
		//std::cout<<"matr a: \n"; print_matrix_spv(a,n,n,5);
		//std::cout<<"matr inv: \n"; print_matrix_spv(inv,n,n,5);
       
	}
	 
	return 0;
    
}


void synchronize(int total_threads) {  //Bgch, page 180
    
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    
    pthread_mutex_lock(&mutex);
    
    threads_in++;
    if (threads_in >= total_threads) {
        threads_out = 0;
        pthread_cond_broadcast(&condvar_in);
    } else
        while (threads_in < total_threads)
            pthread_cond_wait(&condvar_in,&mutex);
    
    threads_out++;
    if (threads_out >= total_threads) {
        threads_in = 0;
        pthread_cond_broadcast(&condvar_out);
    } else
        while (threads_out < total_threads)
            pthread_cond_wait(&condvar_out,&mutex);
    
    pthread_mutex_unlock(&mutex);
}



/*
int effective_method(double* a,double* inv, int n,double* x_k){
	double s;//,sum; 
	double norm_a,norm_x,norm;	
	double sp,sp2;
	double mat_norm=matrix_norm(n,a);

	for(int i=0;i<n-1;++i){
		s=0.0;
		for (int j=i+1;j<n;++j){
			s+=(a[i*n+j]*a[i*n+j]);
		}
		norm_a=sqrt(s+(a[i*n+i]*a[i*n+i]));
		x_k[i]=a[i*n+i]-norm_a;        
		norm_x=sqrt(x_k[i]*x_k[i]+s);
		if (norm_x<mat_norm * 1e-15){a[i*n+i]=norm_a;  continue; }        
		norm=1.0/norm_x;
		x_k[i]*=norm;
		for (int j=i+1;j<n;++j){
			x_k[j]=a[i*n+j]*norm;
		}

		for (int j=0;j<n;++j){
			sp=0.0;
			sp2=0.0;
			for (int k=i;k<n;++k){
				sp+=x_k[k]*a[j*n+k];
				sp2+=x_k[k]*inv[k*n+j];

			}
			for (int k=i;k<n;++k){
				a[j*n+k]-=(2*sp*x_k[k]);
				inv[k*n+j]-=(2*sp2*x_k[k]);
				//sp2+=a[j*n+i]*inv[j*n+k];			
			}
		}

	}
	//std::cout<<"matr a: \n"; print_matrix_spv(a,n,n,4);
	if (is_zero_matrix(n,a)==-1) {std::cout<<"Info:   Singular matrix!"<<std::endl; return -1;}
	//std::cout<<"Before Gauss:\n";
    //std::cout<<"matr a: \n"; print_matrix_spv(a,n,n,4);
    //std::cout<<"matr inv: \n"; print_matrix_spv(inv,n,n,4);	
    
    
	if (gauss_back_run(a,inv,n)==-1) return -1;
       

    //std::cout<<"After Gauss:\n";
    //std::cout<<"matr a: \n"; print_matrix_spv(a,n,n,4);
	//std::cout<<"matr inv: \n"; print_matrix_spv(inv,n,n,4);	

	
	 return 0;
}
*/



int effective_method(double* a,double* inv, int n,double* x_k,int thread_num, int total_threads,int& flag){


	double eps;
	//if (thread_num==0) std::cout<<"now cnt = "<<cnt<<std::endl;
	if (thread_num==0) eps=1e-15*matrix_norm(n,a);
	//cnt--;
	double s;
	//short flag2=1;
	double norm_a,norm_x,norm;	
	double sp,sp2;
	//double mat_norm = matrix_norm(n,a);

	
	
    int first_index, last_index,tmp;
	//if (n%total_threads==0) tmp= n/total_threads;
	//else tmp= n/total_threads+1;
	tmp= n/total_threads;
    first_index=tmp;
	
        
    first_index=tmp*thread_num;
    last_index=(thread_num >= total_threads - 1)?n:(tmp*thread_num+tmp);
	if ((n%total_threads) && (thread_num==total_threads-1)) last_index=n;
	//std::cout<<thread_num<<": FIRST INDEX = "<<first_index<<" AND LAST_INDEX = "<<last_index<<std::endl;

	for(int i=0;i<n-1;++i){
	    synchronize(total_threads);
        if (thread_num==0){
            s=0.0;
            for (int j=i+1;j<n;++j){
                s+=(a[i*n+j]*a[i*n+j]);
            }
            norm_a=sqrt(s+(a[i*n+i]*a[i*n+i]));
            x_k[i]=a[i*n+i]-norm_a;        
            norm_x=sqrt(x_k[i]*x_k[i]+s);
            //if (norm_x<eps) { a[i*n+i]=norm_a; continue; }        
			//if (norm_x>eps){	
			if (norm_x>eps) norm=1.0/norm_x;
			x_k[i]*=norm;
			for (int j=i+1;j<n;++j){
				x_k[j]=a[i*n+j]*norm;
			}
			
        }
        synchronize(total_threads);
        /*
            f1___l1_____________  
            |____|____|____|____| (n)

                f2__l2__________  
                |___|___|___|___| (n)
        */
        
		//std::cout<<"\nTHREAD_NUM ======== "<<thread_num<<std::endl;
		for (int j=first_index;j<last_index;j++){
			//std::cout<<": FIRST INDEX = "<<first_index<<" AND LAST_INDEX = "<<last_index<<" andj = "<<j<<std::endl;
			sp=0.0;
			sp2=0.0;
			for (int k=i;k<n;++k){
				sp+=x_k[k]*a[j*n+k];
				sp2+=x_k[k]*inv[k*n+j];
			}
			for (int k=i;k<n;++k){
				a[j*n+k]-=(2*sp*x_k[k]);
				inv[k*n+j]-=(2*sp2*x_k[k]);	
			}
		
		}
		//synchronize(total_threads);	
		//std::cout<<"\nleave THREAD_NUM = "<<thread_num<<std::endl;
	}
	double tmpp=1;
	synchronize(total_threads);
	//---------------Gauss-------------------------
	//  _______________________
	// |*			||*	  *   * |
	// |	*		||*	  *   * |
	// |		*	||*	  *   * |
	// |____________||__________|


	for (int i=n-1;i>0;--i){
		for (int j=i-1;j>=0;--j){
			tmpp=a[i*n+j];
			//a[i*n+j]=0.0;
			for (int k=first_index;k<last_index;k++){
				inv[j*n+k]-=tmpp*inv[i*n+k]/a[i*n+i];
			}
		}
	}

	synchronize(total_threads);
		//if (thread_num==0){
					//std::cout<<"matr inv: \n"; print_matrix_spv(inv,n,n,4);
					//std::cout<<"matr a: \n"; print_matrix_spv_row(inv,n,n,4);
		//		}  
	
	if (thread_num==0){
		for (int i=0;i<n;++i){
			if (std::abs(a[i*n+i])<eps) flag=0;
			tmpp=a[i*n+i];
			for (int j=0;j<n;++j){
				//a[j*n+i]=a[j*n+i]/tmpp;
				inv[i*n+j]=inv[i*n+j]/tmpp;
			}
		}
	}
	
	return 0;
}


//183
void* effective_method(void* pa){
    ARGS* pargs = (ARGS*)pa;
    //std::cout<<"Info:   thread "<< pargs->thread_num<<"started\n";
   // synchronize(pargs -> total_threads);
    double t=get_full_time();
    effective_method(pargs->matrix,pargs->inv,pargs -> n,pargs -> x_k, pargs -> thread_num, pargs-> total_threads, pargs->flag_error);
    synchronize(pargs -> total_threads);
    
    pargs->time = get_full_time()-t;
    //std::cout<<"\nInfo from function: TIME = "<<pargs->time<<std::endl;

    return 0;
}

void* residual(void *pa){
    ARGS* pargs = (ARGS*)pa;
	double ma = 0.0;
	double str_sum = 0.0;
	double sum=0.0;
	pargs->nor = 0;
    int p = pargs->total_threads;
	double *a = pargs->matrix_orig;
    double *inv = pargs->inv;
    int n = pargs->n;

    mtx.lock();            		//  will not allow 
    (pargs->thread_num)++;		//  other threads 		
    int k = pargs->thread_num;  //  to enter until         
    mtx.unlock();               //  it is completed

    
    for ( int i = k-1; i < n; i+=p){
        str_sum = 0.0;
        for (int j = 0; j < n; j++){
            sum=0.0;
            for (int y = 0; y < n; y++){
                sum += (a[y * n + i] * inv[y * n + j]);
            }
            if (i == j) sum -= 1.0;
            str_sum+=fabs(sum);
        }
        pargs->ArrayForNorm[i] = fabs(sum);
    }
    synchronize(p);
    if (k == 1){
        for (int i = 0; i < n; i++){
            if (pargs->ArrayForNorm[i] > ma){
                ma= pargs->ArrayForNorm[i];
            }
        }
        pargs->nor = ma;
    }
    return 0;
}
