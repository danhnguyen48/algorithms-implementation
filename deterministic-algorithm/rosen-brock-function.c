#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define tol 1.0e-10
#define iters 10000
#define initilizedAlpha 3.0
#define MAX_DOUBLE 100000000
#define small_number 1.0e-16	// detect possible numerical errors caused by small divisions
#define sigma 1.0e-4			// backtracking parameter
#define rho	0.7			        // backtracking parameter
#define lev_mar_parameter 2				// Levenberg-Marquardt ni parameter

double * gradientF(double x1, double x2) {

    double * gradientVector = (double *) malloc(2*sizeof(double));
    *gradientVector = -400 * (x2 - x1*x1) * x1 - 2 * (1 - x1);
    *(gradientVector + 1) = 200 * (x2 - x1*x1);

    return gradientVector;

}

double * hessianMatrix(double x1, double x2) {

    // hessian[0]: dx1^2, hessian[1]: dx2^2, hessian[2]: dx1dx2
    double * hessianMatrix = (double *) malloc(3*sizeof(double));
    *hessianMatrix = 1200*x1*x1 - 400*x2 + 2;
    *(hessianMatrix + 1) = 200;
    *(hessianMatrix + 2) = -400 * x1;

    return hessianMatrix;

}

double fCalculation(double x1, double x2) {

    return 100 * (x2 - x1*x1) * (x2 - x1*x1) + (1 - x1) * (1 - x1);

}

double backtracking(double * x, double * gradient) {
    
    double alpha = initilizedAlpha;
    double val = fCalculation(x[0], x[1]);
    double *J = gradientF(x[0], x[1]);
    double *xx = malloc(2*sizeof(double));
    double val1;

    *xx = x[0] + alpha*(-gradient[0]); *(xx + 1) = x[1] + alpha*(-gradient[1]);
    val1 = fCalculation(xx[0], xx[1]);

    while (val1 > val + sigma*alpha*(J[0]*(-gradient[0]) + J[1]*(-gradient[1]))) {
        alpha = alpha * rho;
        *xx = x[0] + alpha*(-gradient[0]);
        *(xx + 1) = x[1] + alpha*(-gradient[1]);
        val1 = fCalculation(xx[0], xx[1]);
    }

    free(J); free(xx);

    return alpha;

}

int normalize(double * v){

	// normalizes the 2 dimensional vector v. Returns 1 if a division by a very small number is performed
	int div = 0;
	double norm = sqrt(v[0]*v[0] + v[1]*v[1]);
	if (norm < small_number) div = 1;

	v[0] = v[0]/norm;
	v[1] = v[1]/norm;
	return div;

}

void steepest_A_modified_step(double * x, double * g, double * new){

	// Steepest descend and modified newton step
	normalize(g);
	double alpha=backtracking(x, g);

	new[0]=x[0]-alpha*g[0];
	new[1]=x[1]-alpha*g[1];

}

int solve_system(double * A, double * b, double * result){
	// Solve a 2*2 linear system Ax=b, where A is symmetric, using Gauss Elimination. Returns 1 if a div by a small number is performed
	int div=0;
	if (fabs(A[0])<small_number) div=1;

	double pivot=A[2]/A[0];
	double coeff=A[1]-pivot*A[2];
	double grad=b[1]-pivot*b[0];

	result[1]=grad/coeff;
	result[0]=(b[0]-A[2]*result[1])/A[0];

	return div;
}

void stepX(double * x, double * delta, double * new){

	// Original Newton and Levenberg-Marquardt step
	new[0]=x[0]-delta[0];
	new[1]=x[1]-delta[1];

}

double * JtJ(double * x){

	//result will contain the result of (Jacobian^T)*(Jacobian) evaluation at x
    double *result = malloc(3*sizeof(double));
	result[0] = 400 * x[0] * x[0] + 1;
	result[1] = 100;
	result[2] = -200 * x[0];

    return result;

}

double * Jtf(double * x){

	//calculates J^T*f in the L-M algorithm notation
    double *result = malloc(2*sizeof(double));
	result[0] = -(200 * x[0] * (x[1]-x[0]*x[0]) - x[0] + 1);
	result[1] = 100 * (x[1]-x[0]*x[0]);

	return result;
    
}

double max(double* matrix){
	// Returns the bigger diagonal element of a 2*2 matrix
    return matrix[1] > matrix[0] ? matrix[1] : matrix[0];
}

double * damp_sum(double damp, double* matrix){

	// Calculates matrix+damp*(diag(matrix))
    double *sum = malloc(3*sizeof(double));
	sum[0] = matrix[0] + damp*matrix[0];
	sum[1] = matrix[1] + damp*matrix[1];
	sum[2] = matrix[2];

	return sum;

}

int main(int args, char *argv[]) {

    int iter = 0;
    double tolCalculated = MAX_DOUBLE;
    double x1 = -1.5, x2 = -1; // Initialize point

    double * x_k = malloc(2*sizeof(double));
    double * x_k1 = malloc(2*sizeof(double));
    double * gradient = malloc(2*sizeof(double));
    double * hessian = malloc(3*sizeof(double));
    double * d_k = malloc(2*sizeof(double));

    // Steepest Descent
    *x_k = x1; *(x_k + 1) = x2;
    while (iter <= iters && tolCalculated >= tol) {

        gradient = gradientF(x_k[0], x_k[1]);
        steepest_A_modified_step(x_k, gradient, x_k1);

        tolCalculated = fCalculation(x_k[0], x_k[1]) - fCalculation(x_k1[0], x_k1[1]);
        double *temp = x_k; x_k = x_k1; x_k1 = temp;

        iter++;

    }
    printf("STEEPEST: Amount of iters: %d, f: %.16f, x1: %.12f, x2: %.12f\n", iter, fCalculation(x_k[0], x_k[1]), x_k[0], x_k[1]);

    // Original Newton
	iter = 0;       				
	tolCalculated = MAX_DOUBLE;
	*x_k = x1; *(x_k + 1) = x2;
	while (iter <= iters && tolCalculated >= tol) {

		gradient = gradientF(x_k[0], x_k[1]);
        hessian = hessianMatrix(x_k[0], x_k[1]);
		solve_system(hessian, gradient, d_k);
		stepX(x_k, d_k, x_k1);
	
        tolCalculated = gradient[0]*gradient[0] + gradient[1]*gradient[1];
		double *temp = x_k; x_k = x_k1; x_k1=temp;

		iter++;

	}
	printf("NEWTON: Amount of iters: %d, f: %.16f, x1: %.12f, x2: %.12f\n", iter, fCalculation(x_k[0], x_k[1]), x_k[0], x_k[1]);

    // Modified Newton
    iter = 0;
	tolCalculated = MAX_DOUBLE;
	*x_k = x1; *(x_k + 1) = x2;
	while (iter <= iters && tolCalculated >= tol) {
		
        gradient = gradientF(x_k[0], x_k[1]);
        hessian = hessianMatrix(x_k[0], x_k[1]);
		solve_system(hessian, gradient, d_k);
		steepest_A_modified_step(x_k, d_k, x_k1);

		tolCalculated = gradient[0]*gradient[0] + gradient[1]*gradient[1];
		double *temp = x_k; x_k = x_k1; x_k1=temp;

		iter++;

	}
	printf("MODIFIED NEWTON: Amount of iters: %d, f: %.16f, x1: %.12f, x2: %.12f\n", iter, fCalculation(x_k[0], x_k[1]), x_k[0], x_k[1]);

    //LEVENBERG-MARQUARDT
	iter = 0;
	tolCalculated = MAX_DOUBLE;
	*x_k = x1; *(x_k + 1) = x2;
	
	double * delta = malloc(2*sizeof(double));
	double * r_side = malloc(2*sizeof(double));
	double * diag_sum = malloc(3*sizeof(double));
    double * JtJResult = malloc(3*sizeof(double));

	JtJResult = JtJ(x_k);
	double damp=max(JtJResult);
	
	while (iter <= iters && tolCalculated >= tol) {

		JtJResult = JtJ(x_k);
		diag_sum = damp_sum(damp, JtJResult); //calculation of JtJ+damp*(diag(JtJ))
		r_side = Jtf(x_k); //calculation of Jtf
		
		solve_system(diag_sum, r_side, delta); //calculation of the descent direction
        stepX(x_k, delta, x_k1);

		if (fCalculation(x_k1[0], x_k1[1]) < fCalculation(x_k[0], x_k[1])){
			tolCalculated = fCalculation(x_k[0], x_k[1]) - fCalculation(x_k1[0], x_k1[1]);
			double *temp = x_k; x_k = x_k1; x_k1 = temp;
			damp = damp / lev_mar_parameter;
		} else damp = damp * lev_mar_parameter;

		iter++;	
		
	}
    printf("LEVENBERG-MARQUARDT: Amount of iters: %d, f: %.16f, x1: %.12f, x2: %.12f\n", iter, fCalculation(x_k[0], x_k[1]), x_k[0], x_k[1]);

    return 0;

}