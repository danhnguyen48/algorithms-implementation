#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//algorithm parameters
#define tol			 1.0e-10 		//tolerance parameter for stopping the algorithms
#define n_iter 		 1.0e+4			//maximum number of iterations allowed
#define sigma 		 1.0e-4			//backtracking parameter
#define rho			 0.7			//backtracking parameter
#define lev_par 	 2				//Levenberg-Marquardt ni parameter
#define small_number 1.0e-16		//can be used to detect possible numerical errors caused by small divisions

double f (double* x){
	//rosenbrock function calculator
	return 100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])+(1-x[0])*(1-x[0]);
}

int Df(double* x, double* grad){
	//grad will contain the result of the gradient evaluation at x 
	grad[0]=-400*x[0]*(x[1]-x[0]*x[0])-2*(1-x[0]);
	grad[1]=200*(x[1]-x[0]*x[0]);

	return 0;
}

int Hf(double*x, double* hessian){
	//hessian will contain the result of the hessian evaluation at x
	hessian[0]=2-400*(x[1]-3*x[0]*x[0]); 	//dx1^2
	hessian[1]=200; 						//dx2^2
	hessian[2]=-400*x[0]; 					//dx1dx2

	return 0;

}

int JtJ(double* x, double* result){
	//result will contain the result of (Jacobian^T)*(Jacobian) evaluation at x
	result[0]=400*x[0]*x[0]+1;
	result[1]=100;
	result[2]=-200*x[0];

	return 0;
}

int Jtf(double* x, double* result){
	//calculates J^T*f in the L-M algorithm notation
	result[0]=-(200*x[0]*(x[1]-x[0]*x[0])-x[0]+1);
	result[1]=100*(x[1]-x[0]*x[0]);

	return 0;
}

int normalize(double* v){
	//normalizes the 2 dimensional vector v. Returns 1 if a division by a very small number is performed
	int div=0;
	double norm=sqrt(v[0]*v[0]+v[1]*v[1]);
	if (norm < small_number) div=1;

	v[0]=v[0]/norm;
	v[1]=v[1]/norm;
	return div;
}


int solve_system(double* A, double* b, double* result){
	//solve a 2*2 linear system Ax=b, where A is symmetric, using Gauss Elimination. Returns 1 if a div by a small number is performed
	int div=0;
	if (fabs(A[0])<small_number) div=1;

	double pivot=A[2]/A[0];
	double coeff=A[1]-pivot*A[2];
	double grad=b[1]-pivot*b[0];

	result[1]=grad/coeff;
	result[0]=(b[0]-A[2]*result[1])/A[0];

	return div;
}

int multiply(double* matrix, double* vector, double* result){
	//calculates the product between a 2*2 symmetric matrix and a 2*1 vector
	result[0]=matrix[0]*vector[0]+matrix[2]*vector[1];
	result[1]=matrix[2]*vector[0]+matrix[1]*vector[1];

	return 0;
}

double max(double* matrix){
	//returns the bigger diagonal element of a 2*2 matrix
	double result=matrix[0];
	if (matrix[1]>matrix[0]) result=matrix[1];

	return result;
}

int damp_sum(double damp, double* matrix, double* sum){
	//calculates matrix+damp*(diag(matrix))
	sum[0]=matrix[0]+damp*matrix[0];
	sum[1]=matrix[1]+damp*matrix[1];
	sum[2]=matrix[2];

	return 0;
}


double backtracking(double* x, double*g){
	//calculates the parameter alpha for each step of the algorithm, g must be a unit vector pointing to the opposite of the descent direction
	double alpha=3;
	double val=f(x);
	double* J=malloc(2*sizeof(double));
	double* copy=malloc(2*sizeof(double));
	double f1=0;

	Df(x, J);
	copy[0]=x[0]+alpha*(-g[0]); copy[1]=x[1]+alpha*(-g[1]);
	f1=f(copy);
	
	while(f1>val+sigma*alpha*(J[0]*(-g[0])+J[1]*(-g[1]))){
		alpha=alpha*rho;
		copy[0]=x[0]+alpha*(-g[0]); copy[1]=x[1]+alpha*(-g[1]);
		f1=f(copy);
	}

	free(J); free(copy);
	
	return alpha;
}


int mod_step(double* x, double* g, double* new){
	//steepest descend and modified newton step
	normalize(g);
	double alpha=backtracking(x, g);

	new[0]=x[0]-alpha*g[0];
	new[1]=x[1]-alpha*g[1];

	return 0;
}

int step(double* x, double* delta, double* new){
	//classical Newton and Levenberg-Marquardt step
	new[0]=x[0]-delta[0];
	new[1]=x[1]-delta[1];

	return 0;
}

int main(int argc, char *argv[]){

	double x0[2]={-1.5, -1};

	int iter=0;
	double diff=1;


	double* x=malloc(2*sizeof(double));
	double* x_succ=malloc(2*sizeof(double));
	double* temp=NULL; //to keep x and x_succ pointing to different memory locations during execution
	double* grad=malloc(2*sizeof(double));
	double* matrix=malloc(3*sizeof(double));
	double* dir=malloc(2*sizeof(double));

	//STEEPEST DESCEND
    iter=0;       				
	diff=1;
	x[0]=x0[0]; x[1]=x0[1];
	while (iter<n_iter && diff>=tol){
		Df(x, grad);
		mod_step(x, grad, x_succ);  

		diff=f(x)-f(x_succ);	
		temp=x; x=x_succ; x_succ=temp; 
		iter++;
	}
	printf("STEEPEST: n_iter= %d, f= %.12f, x*= (%.12f, %.12f)\n", iter, f(x), x[0], x[1]);

	//CLASSICAL NEWTON
	iter=0;       				
	diff=1;
	x[0]=x0[0]; x[1]=x0[1];
	
	
	while (iter<n_iter && diff>=tol){
		Df(x, grad); 
		Hf(x, matrix); 	//evaluation of the hessian at x
		solve_system(matrix, grad, dir);
		step(x, dir, x_succ);
	
		diff=grad[0]*grad[0]+grad[1]*grad[1];
		temp=x; x=x_succ; x_succ=temp;
		iter++;
	}
	printf("NEWTON: n_iter= %d, f= %.12f, x*= (%.12f, %.12f)\n", iter, f(x), x[0], x[1]);

	//MODIFIED NEWTON-BACKTRACKING
	iter=0;
	diff=1;
	x[0]=x0[0]; x[1]=x0[1];
	
	
	while (iter<n_iter && diff>=tol){
		Df(x, grad); 
		Hf(x, matrix); 	
		solve_system(matrix, grad, dir);
		mod_step(x, dir, x_succ);

		diff=grad[0]*grad[0]+grad[1]*grad[1];
		temp=x; x=x_succ; x_succ=temp;
		iter++;
	}
	printf("MODIFIED NEWTON: n_iter= %d, f= %.12f, x*= (%.12f, %.12f)\n", iter, f(x), x[0], x[1]);

	//LEVENBERG-MARQUARDT
	iter=0;
	diff=1;
	x[0]=x0[0]; x[1]=x0[1];
	
	double* delta=malloc(2*sizeof(double));
	double* r_side=malloc(2*sizeof(double));
	double* diag_sum=malloc(3*sizeof(double));

	JtJ(x, matrix);
	double damp=max(matrix);
	
	while (iter<n_iter && diff>=tol){

		JtJ(x, matrix);		
		damp_sum(damp, matrix, diag_sum); //calculation of JtJ+damp*(diag(JtJ))
		Jtf(x, r_side);  			 //calculation of Jtf
		
		solve_system(diag_sum, r_side, delta);		//calculation of the descent direction
		step(x, delta, x_succ);

		if (f(x_succ)<f(x)){
			diff=f(x)-f(x_succ);
			temp=x; x=x_succ; x_succ=temp;
			damp=damp/lev_par;
		} else damp=damp*lev_par;

		iter++;	
		
	}
	printf("LEVENBERG-MARQUARDT: iter= %d, f= %.12f, x*=(%.12f, %.12f)\n", iter, f(x), x[0], x[1]);
	
	//free the allocated memory
	free(x); free(x_succ); free(grad); free(matrix); free(dir); free(delta); free(r_side); free(diag_sum);

	return 0;
}