#include "func.hpp"


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
	
	if      (k==1) return n-p+1;
	else if (k==2){
        if (i==j) return 2;
        else if ((p-q)==1) return -1;
        else return 0;

    }
	else if (k==3){
        if ((i==j) && (j<n)) return 1;
        else if (j==n) return i;
        else if (i==n) return j;
        else return 0;
    }
	else if (k==4) return 1.0/(p+q-1);
	return 0;
}

double matrix_norm(int n, double *a)
{
	double sum=0;
    double max=-10000;
	double norm_m;

   	for(int i = 0; i < n; i++)
    { 
        for(int j = 0; j < n; j++)
            sum += std::abs(a[j*n+i]); 
        if(max < sum)
            max = sum;
        sum = 0.0; 
    } 
    norm_m = max; 
    return norm_m;
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
				return -2;
			}
			
			else{
				while(!feof(f)){
					q=fscanf(f,"%s",&r);
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
							return -2;
						}
						arr[j*n+i]=a;
					}
				}
				fclose(f);
			}
			double norm = matrix_norm(n,arr);
			//if(is_zero_matrix(n, arr)==-1) return -1;
			for (int i=0;i<n*n;++i){ if (std::abs(arr[i])<=1e-15*norm){ cnt++;}  }
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

int sym_matr(double* a, int n){
	double norm = matrix_norm(n, a);
	for(int i = 0; i < n; ++i)
	{
		for(int j = i; j < n; ++j)
		{
			if(std::abs(a[i * n + j] - a[j * n + i]) > norm*1e-15) { return -1; }
		}
	}
	return 1;
}

int to_tridiag_form(double* a, double* y,double* x_k, double* z,int n){ 
    double s;
    double norm_a,norm_x,norm;
    double alpha=0.0;
	double mat_norm=matrix_norm(n,a);
	

    for(int i=0;i<n-2;++i){
	
		s=0.0;
		for (int j=i+2;j<n;++j){
			s+=(a[i*n+j]*a[i*n+j]);
		}
		norm_a=sqrt(s+(a[i*n+i+1]*a[i*n+i+1]));

        //-----------x_k----------
		if (s<1e-15*mat_norm){                                        //подумать тут ещё 
			a[i * n + (i + 1)] = norm_a;
	    	a[(i + 1) * n + i] = norm_a;

			for(int k = i + 2; k < n; k++) {
		    	a[i * n + k] = 0.0;
		    	a[k* n + i] = 0.0;
	    	}
			continue;
		}

		x_k[i+1]=a[i*n+i+1]-norm_a;  
		norm_x=sqrt(x_k[i+1]*x_k[i+1]+s);
		if (norm_x<1e-15*mat_norm){ 
			continue;
		}		
		norm=1.0/norm_x;
		x_k[i+1]*=norm;

        //std::cout<<"x = (";
        //std::cout<<x_k[i+1]<<" ";
		for (int j=i+2;j<n;++j){
			x_k[j]=a[i*n+j]*norm;
            //std::cout<<x_k[j]<<" ";
		}
        //std::cout<<")\n";    

        //------------y-----------
        //std::cout<<"y = (";
        for (int k=i+1;k<n;++k){
            y[k]=0.0;
            for (int j=i+1;j<n;++j){
                y[k]+=a[j*n+k]*x_k[j];
            }
            //std::cout<<y[k]<<" ";
        }
        //std::cout<<")\n";
        
        //-----alpha = 2(x,y)-----
		alpha=0.0;
        for (int k=i+1;k<n;++k){
            alpha+=x_k[k]*y[k];
        }
        alpha*=2.0;
        
        //-----z = 2y-alpha*x-----
        //std::cout<<"z = (";
        for (int k=i+1;k<n;++k){
            z[k]=2.0*y[k]-alpha*x_k[k];
            //std::cout<<z[k]<<" ";
        }
        //std::cout<<")\n";

        //---matrB = A-zx*-xz*----

        for (int k=i+1;k<n;++k){
            for (int j=i+1;j<n;++j){
                a[j*n+k]-=z[k]*x_k[j]+x_k[k]*z[j];   
            }
        }

		a[i * n + (i + 1)] = norm_a;
	    a[(i + 1) * n + i] = norm_a;

	    for(int k = i + 2; k < n; k++) {
		    a[i * n + k] = 0.0;
		    a[k* n + i] = 0.0;
	    }
        //std::cout<<"MatrB: \n\n"; print_matrix_spv(a,n,n, n);

    }
	if (std::abs(a[(n-1)*n+n-2])<1e-15*mat_norm) a[(n-1)*n+n-2] = 0.0;   
	if (std::abs(a[(n-2)*n+n-1])<1e-15*mat_norm) a[(n-2)*n+n-1] = 0.0;             
    return 0;

}

double search_m_eps(){
	double eps=1;
	while (1+eps/2>1){
		eps=eps/2.0;
	}

	return eps;
}

int sign_changes(double* a,int n,double lmbd){
	int m=0;
	double alpha;
	double max_a=-128, max_b=-128;
	int tmp;
	double tmp2,tmp3;
	double x,y,u,v,g;
	double m_eps=search_m_eps();
	double a1,b1;


	//---------------alpha=4*max{ max{|a_i|}, max{|b_i|} }----------------
	for (int i=0;i<n-1;++i){
		tmp=i*n+i;
		tmp2=std::abs(a[tmp]-lmbd);
		tmp3=std::abs(a[tmp+1]);
		if (tmp3<m_eps) a[tmp+1]=0.0;
		if (tmp2>max_a) max_a=tmp2;
		if (tmp3>max_b) max_b=tmp3; 
	}
	//(n-2)*n+(n-2)=n*n-2n+n-2=n*n-n-2
	//(n-1)*n+(n-1)=n*n-n+n-1=n*n-1=n*n-n-2+n+1
	if (std::abs(a[tmp+n+1]-lmbd)>max_a) max_a=std::abs(a[tmp]-lmbd);

	alpha = (max_a>max_b)?max_a:max_b;
	alpha*=4.0;
	double aalpha=1.0/alpha;

	//for (int i=0;i<n*n;++i) a[i]*=aalpha;

	x=aalpha*(a[0]-lmbd);
	y=1.0;
	if (x<0) m++;

	double m_eeps=1/m_eps;
	double ma,mma;

	for (int k=1;k<n;++k){
		
		tmp=k*n+k;
		a1=aalpha*(a[tmp]-lmbd);
		b1=aalpha*(a[tmp-1]);
	
		tmp2=std::abs(x);
		tmp3=std::abs(b1*b1*y);
		
		ma=(tmp2>tmp3)?tmp2:tmp3;
		mma=1.0/ma;

		g=m_eeps*mma;

		u=g*(a1*x-b1*b1*y);
		v=g*x;

		if (u*x<0.0) m++;

		x=u;
		y=v;


	}
	return m;
}

int search_values(int n,double* matr,double* lmbd_values,double eps){
	double a0, b0;
	double a,b,c;
	double m_eps=search_m_eps();
	double mat_norm=matrix_norm(n,matr);
	b0=mat_norm+m_eps;
	a0=b0*(-1.0);
	//std::cout<<"a0 = "<<a0<<" & b0 = "<<b0<<std::endl;
	//lmbd_values[0]=a0;
	int diff=0;
	a=a0;
	b=b0;
	int k=0;
	while (k<n){
		
		std::cout<<"Iteration "<<k<<std::endl;
		while (b-a>eps){
			c=(a+b)*0.5;
			//std::cout<<"a = "<<a<<std::endl;
			//std::cout<<"c = "<<c<<std::endl;
			//std::cout<<"b = "<<b<<std::endl;
			if (sign_changes(matr,n,c)<k+1) a=c;
			else b=c;
		}
		c=(a+b)*0.5;
		diff=sign_changes(matr,n,b)-sign_changes(matr,n,a);
		if (diff<0){std::cout<<"o la la "<<std::endl; return -1;}
		std::cout<<"diff = "<<diff<<" and c = "<<c<<std::endl;
		for (int i=0;i<diff;++i) {
			lmbd_values[k+i]=c;
			if (std::abs(c)<m_eps*mat_norm) lmbd_values[k+i]=0.0;// подумать тут ещё...
		}
		k+=diff;
		a=c;
		b=b0;		
	}

	
	//std::cout<<"n_(a) = "<<sign_changes(matr,n,a)<<std::endl;
	//std::cout<<"n_(b) = "<<sign_changes(matr,n,b)<<std::endl;

	return 1;
}


int discrepancy_inv1(int n,double* a,double* lmbd_values){
	double matr_norm;
	matr_norm = matrix_norm(n,a);
	double s = 0;

	//i*n+i
	//(i+1)*n+i+1=i*n+n+i=(i*n+i)+n+1 
	for(int i = 0; i < n; i ++){
		s += (a[i * n + i] - lmbd_values[i]);
	}
	return std::abs(s);
}

void foo(char* a, char* b){
	int i=0;
	while (1){
		if (b[i]==0) break;
		a[i]=b[i];
		i++;
	}

}
