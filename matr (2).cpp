#include "matr.hpp"


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




double norm(int n, double* ar, double* inv) {
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

double matrix_norm(int n, double *a)
{
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








void ineffective_method(double* a,double* inv, int n){
	double t,t_,u;
	int h=0;
	int k;
	int y=0;
	//int nn=n-1;

	double eps= 1.1e-16;
	for (k=0;k<n;++k){
		if (abs(a[k*n+h])<eps){
			y=k;
			while (y<n && abs(a[y*n+h])<eps) y++;
			if (y==n) h++;
			else {
				for (int m=0;m<n;++m){
					t=a[y*n+m];
					a[y*n+m]=a[k*n+m];
					a[k*n+m]=t;
					t_=inv[y*n+m];
					inv[y*n+m]=inv[k*n+m];
					inv[k*n+m]=t_;	
				}
			}
			#ifdef DEBUG
				printf("after changing strings: \n");
				print_matrix(a,n);
				printf("after changing strings unit matrix: \n");
				print_matrix(inv,n);
			#endif
		}

		for (int i=k+1;i<n;++i){
			u=a[i*n+h];
			for (int j=0;j<n;++j){
				#ifdef DEBUG
					printf("a[i][j] = a[i][j]-(a[i-k][h]*a[k][j])/a[k][h];\n");
					printf("%lf = %lf - ( %lf * %lf ) / %lf\n",a[i*n+j],a[i*n+j],u,a[k*n+j],a[k*n+h]);
				
					printf("inv[i][j] = inv[i][j]-(u_*a[k][j])/a[k][h];\n");
					printf("%lf = %lf - ( %lf * %lf ) / %lf\n",inv[i*n+j],inv[i*n+j],u,inv[k*n+j],a[k*n+h]);
				#endif			
				a[i*n+j]=a[i*n+j]-(u*a[k*n+j])/a[k*n+h];
				inv[i*n+j]=inv[i*n+j]-(u*inv[k*n+j])/a[k*n+h];
				#ifdef DEBUG
					printf("after changing: \n");
					print_matrix(a,n);

					printf("after changing unit: \n");
					print_matrix(inv,n);
				#endif
			}
		
		}
		h++;
	}
	//if (check(a,n)==-1) {printf(" Oh no! Singular matrix!\n "); return; };

	for (int i=0;i<n;++i){
		u=a[i*n+i];
		
		for (int j=0;j<n;++j){
			a[i*n+j]=a[i*n+j]/u;
			inv[i*n+j]=inv[i*n+j]/u;

		}
	}
	#ifdef DEBUG
		printf("after changing: \n");
		print_matrix(a,n);
		printf("after changing unit: \n");
		print_matrix(inv,n);
	#endif

	h=n-1;
	for (k=n-1;k>0;--k){
		for (int i=k-1;i>-1;--i){
			u=a[i*n+h];
			for (int j=n-1;j>-1;--j){
				#ifdef DEBUG
					printf("a[i][j] = a[i][j]-(a[i-k][h]*a[k][j])/a[k][h];\n");
					printf("%lf = %lf - ( %lf * %lf ) / %lf\n",a[i*n+j],a[i*n+j],u,a[k*n+j],a[k*n+h]);
				
					printf("inv[i][j] = inv[i][j]-(u_*a[k][j])/a[k][h];\n");
					printf("%lf = %lf - ( %lf * %lf ) / %lf\n",inv[i*n+j],inv[i*n+j],u,inv[k*n+j],a[k*n+h]);
				#endif			
				a[i*n+j]=a[i*n+j]-(u*a[k*n+j])/a[k*n+h];
				inv[i*n+j]=inv[i*n+j]-(u*inv[k*n+j])/a[k*n+h];
				
				#ifdef DEBUG
					printf("after changing: \n");
					print_matrix(a,n);

					printf("after changing unit: \n");
					print_matrix(inv,n);
				#endif
			}
		
		}
		h--;
	}

	//printf("inverse matrix: \n");
	//print_matrix(inv,n);


}
