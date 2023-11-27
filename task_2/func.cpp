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

int SymMatr(double* a, int n){
	double norm = matrix_norm(n, a);
	for(int i = 0; i < n; ++i)
	{
		for(int j = i; j < n; ++j)
		{
			if(std::abs(a[i * n + j] - a[j * n + i]) > norm*1e-15) {std::cout<<"Info:   the matrix is not symmetrical!\n"; return -1; }
		}
	}
	return 1;
}

int to_tridiag_form(double* a, double* y,double* x_k, double* z,double* matrB,int n){ //досчитать!*************
    double s;
    double norm_a,norm_x,norm;
    double alpha=0.0;
    
    //std::cout<<"MatrB: \n\n"; print_matrix_spv(matrB,n,n, n);
    for(int i=0;i<n-2;++i){
        //std::cout<<"ITERATION "<<i<<std::endl;
		s=0.0;
		for (int j=i+2;j<n;++j){
			s+=(a[i*n+j]*a[i*n+j]);
		}
        //std::cout<<"SUM! = "<<s<<std::endl;
		norm_a=sqrt(s+(a[i*n+i+1]*a[i*n+i+1]));
        //std::cout<<"norm_a = "<<norm_a<<std::endl;
        //-----------x_k----------
		x_k[i+1]=a[i*n+i+1]-norm_a;      
        //std::cout<<"x_k[i+1] = "<<x_k[i+1]<<std::endl;  
		norm_x=sqrt(x_k[i+1]*x_k[i+1]+s);
	    if (norm_x< /*mat_norm * */ 1e-15){	a[i*n+i+1]=norm_a;  continue; }       //подумать тут ещё 
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
                //std::cout<<"y["<<k<<"]+=a["<<j*n+k<<"]*x_k["<<j<<"] = "<<y[k]<<" += "<<a[j*n+k]<<"*"<<x_k[j]<<std::endl;
            }
            //std::cout<<y[k]<<" ";
        }
        //std::cout<<")\n";
        
        //-----alpha = 2(x,y)-----
        for (int k=i+1;k<n;++k){
            alpha+=x_k[k]*y[k];
        }
        alpha*=2.0;
        //std::cout<<"ALPHA = "<<alpha<<std::endl;
        
        //-----z = 2y-alpha*x-----
        //std::cout<<"z = (";
        for (int k=i+1;k<n;++k){
            z[k]=2.0*y[k]-alpha*x_k[k];
            //std::cout<<z[k]<<" ";
        }
        //std::cout<<")\n";

        //---matrB = A-zx*-xz*----
        
        //matrB[i*n+i]=a[i*n+i];
        
        matrB[i * n + (i + 1)] = norm_a;
	    matrB[(i + 1) * n + i] = norm_a;

	    for(int k = i + 2; k < n; k++) {
		    matrB[i * n + k] = 0.0;
		    matrB[k* n + i] = 0.0;
	    }
        


        //std::cout<<"MatrB: \n\n"; print_matrix_spv(matrB,n,n, n);
        //std::cout<<"MatrB: \n\n "; print_matrix_spv(matrB,n,n, n);

        for (int k=i+1;k<n;++k){
            for (int j=i+1;j<n;++j){
                //std::cout<<matrB[k*n+j]<<" -= "<<z[k]<<" * "<<x_k[j]<<"+"<<x_k[k]<<"*"<<z[j]<<std::endl;
                matrB[k*n+j]-=(z[k]*x_k[j]+x_k[k]*z[j]);   
            }
        }
        std::cout<<"MatrB: \n\n"; print_matrix_spv(matrB,n,n, n);

    }
    return 0;

}

double SearchEps(){
	double eps=1;
	while (1+eps/2>1){
		eps=eps/2.0;
	}

	return eps;
}

int n_(double * matr, double lmbd, int n){
	int cnt = 0;
    double a,b;
    double gamma,alpha;
    double u,v;
    double x,y;
	double m_a = -1, m_b = -1;
    double m_eps=SearchEps();
    lmbd=0.0;// точно нужно здесь?, потом уберу мб..


    //------------------alpha=4max{max(|a_i|),max(|b_i|)}-------------
	for (int i = 0; i < n; ++i){
		if (std::abs(matr[i * n + i]-lmbd) > m_a) 
			m_a = std::abs(matr[i * n + i] - lmbd);
	}
    for (int i = 0; i < n-1; i ++){
		if (std::abs(matr[i * n + i+1]-lmbd) > m_b) 
			m_b = std::abs(matr[i * n + i+1] - lmbd);
    }
	alpha=(m_a>m_b)?m_a:m_b;
	alpha*=4.0; 
        

	if(alpha > 0)
		alpha = 1.0 / alpha;
	else{
		std::cout << "is it possible?" << std::endl;
		return -1;
	}
	
    matr[0]*=alpha;
    matr[1]*=alpha;
    for (int i=1;i<n-1;++i){
        for (int j=-1;j<2;++j){
            matr[i*n+j+i]*=alpha;
        }
        if (std::abs(matr[i*n+i+1])<1e-16){
            matr[i*n+i+1]=0.0;
            matr[(i+1)*n+i]=0.0;
        } 
    }
    matr[(n-1)*n+n-2]*=alpha;
    matr[(n-1)*n+n-1]*=alpha;
    std::cout<<"MatrBBBBBB: \n\n"; print_matrix_spv(matr,n,n, n);

    x=matr[0];
    y=1;

	if(x * y < 0) 
		cnt++;
    double max;
	for(int k = 1; k < n; k ++)
	{
		a = matr[k*n+k];
		b = matr[k*n+k-1];
		max = (std::abs(x) > std::abs(b * b * y)) ? std::abs(x):std::abs(b*b*y);

		gamma = m_eps / max;
		u = gamma * (a*x-b*b*y);
		v = gamma * x;

		if((u*x)<0) cnt ++;
		x = u;
		y = v;
	}
	return cnt;

}

int recursion(double* matr, int n, double eps,double left, double right, double* lmbd_values,int index){
	int res;
    double c;
    /*
	if(if_rec != 1 || iter > maxit)
		return 1;
	iter++;
    */
	int n_a = n_(matr, left, n);
	int n_b = n_(matr, right, n);

	if(n_b - n_a > 0)
	{
		if (right - left > eps) {
			//if_rec = 1;
			res = recursion(matr, n,eps, left, 0.5 * (left + right), lmbd_values,index); //, index, eps, if_rec, maxit, iter, meps, lam, eigen_krat, global_eigen_num);
			if(res == -1) return -1;
			//if_rec = 1;
			res = recursion(matr, n,eps, 0.5 * (left + right), right, lmbd_values,index); //, index, eps, if_rec, maxit, iter, meps, lam, eigen_krat, global_eigen_num);
			if(res == -1) return -1;
		}
		else{
			c = 0.5 * (left + right);
			for(int i = 0; i < n_b - n_a; ++i) lmbd_values[index + i] = c;
            index += n_b - n_a;
			//lam[*global_eigen_num] = c;
			//eigen_krat[*global_eigen_num] = n_b - n_a;
			//(*global_eigen_num)++;

			 //Добавляем кратность

			if(index >= n)
			{
                std::cout<<"thats all"<<std::endl;
				return 1;
			}
		}
		return 1;
	}
	else
	{
		//if_rec = -1;
        std::cout<<"There is nothing here.."<<std::endl;
		return 1;
	}
    return 0;
}


