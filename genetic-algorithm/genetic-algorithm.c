#include <stdio.h>
#include <stdlib.h>
#include "RKF78.h"
#include "rand1.h"
#include <math.h>

#define n_par 11 // time steps
#define d_par 10 // number of treatments in cure
#define m_par 4 // organs cannot be damaged
#define population_size 30000
#define MAX_DOUBLE 100000000
#define CMAXj 15		
#define CCUMj 127
#define N0 20000
#define lambda 0.336
#define lambda_log_theta 9.284023094951992774680543277273261947643
#define NMax 9.5e+11 // max amount of tumor cells
#define tournament_size 10 // number of selections with replacement in each choice of parents
#define probability_cross 0.2
#define probability_multiple 0.1

static float k_j[d_par] = {0.12, 0.0502, 0.0637, 0.1347, 0.0902, 0.0546, 0.0767, 0.1121, 0.0971, 0.0403};
static float t_i[n_par] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 33.0};
static float eta_kj[m_par][d_par] = 
/*{ //curative parameters 
	{0.0036, 0.0098, 0.0061, 0.0009, 0.0003, 0.0108, 0.0045, 0.0021, 0.0096, 0.0125},
	{0.0063, 0.0082, 0.0062, 0.0062, 0.0083, 0.013, 0.0039, 0.0019, 0.0015, 0.005},
	{0.0129, 0.0018, 0.0116, 0.0021, 0.009, 0.0129, 0.0054, 0.0049, 0.0093, 0.0066},
	{0.0053, 0.0086, 0.0067, 0.0029, 0.0089, 0.0054, 0.0042, 0.0095, 0.0112, 0.0092}
};*/
//palliative parameters
{
	{0.00612, 0.01666, 0.01037, 0.00153, 0.00051, 0.01836, 0.00765, 0.00357, 0.01632, 0.02125},
	{0.01071, 0.01394, 0.01054, 0.01054, 0.01411, 0.0221, 0.00663, 0.00323, 0.00255, 0.0085},
	{0.02193, 0.00306, 0.01972, 0.00357, 0.0153, 0.02193, 0.00918, 0.00833, 0.01581, 0.01122},
	{0.00901, 0.01462, 0.01139, 0.00493, 0.01513, 0.00918, 0.00714, 0.01615, 0.01904, 0.01564}
};

typedef struct {
    double drift_i;
} ODE_Parameters;

void gompertz(double t, double N, double *der, void *Params){
	*der = ((N < 1.e-16) ? 0.0 : N*(lambda_log_theta - lambda*log(N) - ((ODE_Parameters *) Params)->drift_i));
}

void mod_gompertz(double t, double N, double *der, void* Params){
	*der = ((N < 1.e-16) ? 0.0 : N*(lambda_log_theta - lambda*log(N)));
}

void init_population(unsigned char* population); // Initialize population
unsigned char areContraint23Verified(unsigned char *Cij); // See if Contraint 2 and 3 are verified?
int mutate(unsigned char* child);
void crossover(unsigned char* parent_a, unsigned char* parent_b, unsigned char* child1, unsigned char* child2);
void select_with_replacement (unsigned char* pop, double* fitness, unsigned char* parent);
void next_generation(unsigned char *pop1, unsigned char *pop2, double *fitness_value);
double curative_fitness(unsigned char *Cij); // Calculate fitness value
int look_for_curative(unsigned char *pop1, unsigned char *pop2, double *best_fitness, unsigned char *best_individual, double *fitness_values); // Look for curative treatment
double palliative_fitness(unsigned char* Cij); // Calculate palliative fitness value
int look_for_palliative(unsigned char* pop1, unsigned char* pop2, double* best_fitness, unsigned char* best_individual, double* fitness_values);// Look for palliative treatment


int main(int argc, char *argv[]) {

    unsigned char *pop_1=(unsigned char*) malloc(population_size*d_par*(n_par-1)*sizeof(unsigned char));
	unsigned char *pop_2=(unsigned char*) malloc(population_size*d_par*(n_par-1)*sizeof(unsigned char));
	double *fitness_values=(double*) malloc(population_size*sizeof(double));			//contains the ordered fitness values of a generation
	
	unsigned char *bestIndividual = (unsigned char*) malloc(d_par*(n_par-1)*sizeof(unsigned char));	//best individual founded by the algorithm
	double best_fitness = MAX_DOUBLE;
	
	init_population(pop_1);	//initialize the population

	int curative = look_for_curative(pop_1, pop_2, &best_fitness, bestIndividual, fitness_values);

    if (curative==0) { // If found a curative treatment

		printf("A curative treatment has been found. Its doses are \n");
		for (int i=1; i<=d_par*(n_par-1); i++) {
			printf("%d ", bestIndividual[i-1]);
			if (i%d_par==0) printf("\n");
		}
		printf("Its fitness value is %f ", best_fitness);
		printf("\n");

	} else { // If not, must look for a palliative solution

		best_fitness = MAX_DOUBLE;

		int palliative = look_for_palliative(pop_1, pop_2, &best_fitness, bestIndividual, fitness_values);

		if (palliative == 0){

			printf("A palliative treatment has been found. Its doses are \n");
			for (int i=1; i<=d_par*(n_par-1); i++) {
				printf("%d ", bestIndividual[i-1]);
				if (i%d_par==0) printf("\n");
			}
			printf("The lifespan of the patient is %f ", 1/best_fitness);
			printf("\n");

		} else printf("Neither a curative nor a palliative treatment have been found");
	}

}

void init_population(unsigned char* population) {
	//fills the memory pointed by pop with individuals 
	unsigned char *drive = population;
	long int seed = -rand();
	int generated = 0;
	while (drive != population+population_size*d_par*(n_par-1)){ // keep running until drive to the last

		//generate random individual
		for (int i=0; i<d_par*(n_par-1); i++){
			generated = (int) (ran1(&seed)*1.0e+4);
			*(drive + i) = generated % CMAXj;
		}

		drive += d_par*(n_par-1);
	}
}

unsigned char areContraint23Verified(unsigned char *Cij) {
	register unsigned char i, j, k;
	unsigned char npar = n_par - 1, *Cijofi;

	//constraint 2
	for (i=0; i < d_par; i++){ 
		unsigned int ccumj = 0U;
		for (j=0; j < npar; j++) ccumj += *(Cij + j*d_par + i);
		if (ccumj > CCUMj) return 0U;
	}

	//constraint 3
	for (i=0, Cijofi=Cij; i < npar; i++, Cijofi += d_par) {
		for (k=0; k < m_par; k++){ 
			double Cseffk = 0.0;
			for (j=0; j < d_par; j++) Cseffk += eta_kj[k][j] * Cijofi[j];
			if (Cseffk > 1.0) return 0U;
		}
	}
	return 1U;
}

int mutate(unsigned char* child) {
	//bitflip mutation
	long int seed = -rand();
	double p=0;
	for (int i=0; i<d_par*(n_par-1); i++) {
		p = ran1(&seed);
		if (p <= probability_multiple) child[i] = CMAXj - child[i];
	}
	return 1;
}

void crossover(unsigned char* parent_a, unsigned char* parent_b, unsigned char* child1, unsigned char* child2) {
	//uniform crossover between the parents generates the children
	long int seed = -rand();
	double p = 0;

	for (int i=0; i<d_par*(n_par-1); i++){
		p = ran1(&seed);
		
		if (p <= probability_cross){ //no swapping 
			child1[i]=parent_a[i];
			child2[i]=parent_b[i];
		} else { //swapping
			child1[i]=parent_b[i];
			child2[i]=parent_a[i];
			
		}
	}

}

void select_with_replacement (unsigned char* pop, double* fitness, unsigned char* parent) {
	
	long int seed = -rand();
	int selection = ((int) (ran1(&seed)*1.0e+4)) % population_size;
	unsigned char* best = pop + selection*d_par*(n_par-1);
	double fit_best = fitness[selection];
	
	unsigned char* next = NULL;
	for (int i=2; i<=tournament_size; i++){
		//new individual
		selection = ((int)(ran1(&seed)*1.0e+4)) % population_size;
		next = pop + selection*d_par*(n_par-1);
		//if its fitness is better than the best one it becomes the best
		if (*(fitness + selection) < fit_best) {
			fit_best = *(fitness + selection); 
			best = next;
        }
	}

	//copying the selected individual
	for (unsigned char i=0; i<d_par*(n_par-1); i++){
		*(parent+i) = *(best+i);
	}

}

void next_generation(unsigned char *pop1, unsigned char *pop2, double *fitness_value) {
    unsigned char *parent_a = malloc(d_par*(n_par-1)*sizeof(unsigned char));
    unsigned char *parent_b = malloc(d_par*(n_par-1)*sizeof(unsigned char));
    unsigned char *child1 = malloc(d_par*(n_par-1)*sizeof(unsigned char));
    unsigned char *child2 = malloc(d_par*(n_par-1)*sizeof(unsigned char));
    for (int i=0; i<population_size/2; i++){
        
        //parents selection
        select_with_replacement(pop1, fitness_value, parent_a);
        select_with_replacement(pop1, fitness_value, parent_b);

        //children generation
        crossover(parent_a, parent_b, child1, child2);

        //children mutation
        mutate(child1);
        mutate(child2);

        //adding the children to the new generation
        for (unsigned char j=0; j<d_par*(n_par-1); j++){
            pop2[i*d_par*(n_par-1)+j] = child1[j];
            pop2[population_size/2+i*d_par*(n_par-1)+j] = child2[j];
        }	
    }
    //free memory 
    free(parent_a); free(parent_b); free(child1); free(child2);
}

double curative_fitness(unsigned char *Cij) {
	register unsigned char i, j;
	ODE_Parameters GompertzParams;
	double N = N0, t = t_i[0];
	double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
	unsigned char curativecounter = 0U, npar = n_par - 1;
	double integral = 0.0, lastt = t, lastN = N;

	if (!areContraint23Verified(Cij))
        return MAX_DOUBLE;

	

	for (i=0; i < npar; i++){ 
		double tfin = t_i[i+1]; 
		GompertzParams.drift_i = 0.0;
		for (j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

		while (t+h < tfin) {
			RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
			if(N > NMax) return MAX_DOUBLE;
			integral += (lastN + N)*(t - lastt);
			lastt = t; lastN = N;
		}

		do { 
            h = tfin - t;
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
            if(N > NMax) return MAX_DOUBLE;
            integral += (lastN + N)*(t - lastt);
            lastt = t; lastN = N;
		} while (t < tfin);

		if (N < 1000) { 
			curativecounter++; 
			if(curativecounter > 2) return integral/2.0; 
		} else curativecounter = 0U;
	}

	return MAX_DOUBLE;
}

int look_for_curative(unsigned char *pop1, unsigned char *pop2, double *best_fitness, unsigned char *best_individual, double *fitness_values) {

	//control variables
	int n_infeasible_gen = 0;			// number of consecutive infeasible generations
	int n_infeasible_ind = 0;			// number of infeasible individuals in a generation
	int n_non_updating_gen = 0; 	    // number of consecutive generations in which the best individual does not change
	int update_best;			        // Boolean variable, indicates that whether the best individual was changed
	int best_pos = 0;				    // position of the best element of the population in the vector

	int count = 0;
	do {
		printf("Step: %d\n", count);
		n_infeasible_ind = 0;
		update_best = 0; 
		
		//initializing the fitness values, computing the best position, counting the infeasibles
		for(int i=0; i<population_size; i++){
			fitness_values[i] = curative_fitness(pop1 + i*d_par*(n_par-1));
			
			if (fitness_values[i] == MAX_DOUBLE) n_infeasible_ind++; // If it's received max double, it means it's infeasible individual
			else if (fitness_values[i] < *best_fitness){ // Update the best individual

                best_pos = i;
                *best_fitness = fitness_values[i];
                update_best = 1;

            }
		}

		//saving the best individual if it has changed in the generation
		if (update_best == 1) {
            for (unsigned char j=0; j<d_par*(n_par-1); j++) { 
                *(best_individual + j) = *(pop1 + best_pos*d_par*(n_par-1)+j);
            }
        }
			
		//if the generation is infeasible increase n_infeasible_gen, otherwise set it to 0 and check if the best is updated
		if (n_infeasible_ind == population_size) n_infeasible_gen++; 
        else {
			n_infeasible_gen = 0;
			//if the best is updated and the generation is not infeasible the number of non-updating gen is set to 0, otherwise it increases
			if(update_best == 1) n_non_updating_gen = 0; 
            else n_non_updating_gen++;
		}

		printf("number of infeasible individuals: %d\n", n_infeasible_ind);
		printf("n_non_updating: %d\n", n_non_updating_gen);
		printf("number of consecutive infeasible generations: %d\n", n_infeasible_gen);

		//creating the next generation
		next_generation(pop1, pop2, fitness_values);

		unsigned char* temp=pop1; pop1=pop2; pop2=temp;

		count++;
	} while (n_infeasible_gen < 10 && n_non_updating_gen < 4 && count < 100);

	if (n_infeasible_gen >= 10) return 1;

	return 0;
}

double palliative_fitness(unsigned char* Cij) {
	
	// Return 1/lifespan. If the patient is dead, return max double
	register unsigned char i, j;
	ODE_Parameters GompertzParams;
	double N = N0, t = t_i[0];
	double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
	unsigned char npar = n_par - 1;
	double integral = 0.0, lastt = t; 

	if (!areContraint23Verified(Cij)) return MAX_DOUBLE;

	//see if the patient survives until the end of the treatment
	for (i=0; i < npar; i++){ 
		double tfin = t_i[i+1];
		GompertzParams.drift_i = 0.0;
		for (j=0; j < d_par; j++)
			GompertzParams.drift_i += k_j[j] * *(Cij++);

		while (t+h < tfin) {

			RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
			if(N > NMax) return MAX_DOUBLE;
			lastt = t;

		}

		do { 

			h = tfin - t;
			RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
			if (N > NMax) return MAX_DOUBLE;
			lastt = t; 

		} while (t < tfin);

	}

	//if the patient survives calculate the effective lifespan
	while (N < NMax) {
		RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, mod_gompertz);
		lastt = t;
	}

	return 1/lastt;
}

int look_for_palliative(unsigned char* pop1, unsigned char* pop2, double* best_fitness, unsigned char* best_individual, double* fitness_values) {
	//LOOK FOR A PALLIATIVE SOLUTION
	int n_infeasible_gen = 0;			// number of consecutive infeasible generations
	int n_infeasible_ind = 0;			// number of infeasible individuals in a generation
	int n_non_updating_gen = 0; 	    // number of consecutive generations in which the best individual does not change
	int update_best;			        // Boolean variable, indicates that whether the best individual was changed
	int best_pos = 0;				    // position of the best element of the population in the vector

	int count = 0;
	do {
		printf("Step %d\n", count);
		n_infeasible_ind = 0;
		update_best = 0; // false

		//initializing the fitness values, computing the best position, counting the infeasibles
		for (int i=0; i<population_size; i++){

			fitness_values[i] = palliative_fitness(pop1+i*d_par*(n_par-1)); 

			if (fitness_values[i] == MAX_DOUBLE) n_infeasible_ind++; // The patient is dead in the treatment 
			else if (fitness_values[i] < *best_fitness){

				best_pos = i;
				*best_fitness = fitness_values[i];
				update_best = 1;
				
			}

		}

		// Saving the best individual if it has changed in the generation
		if(update_best == 1) {
			for (int j=0; j<d_par*(n_par-1); j++) 
				*(best_individual + j) = *(pop1 + best_pos*d_par*(n_par-1)+j);
		}
			
		// If the generation is infeasible increase n_infeasible_gen, otherwise set it to 0 and check if the best is updated
		// If the best is updated, set n_non_updating_gen to 0, otherwise increase it
		if (n_infeasible_ind == population_size) n_infeasible_gen++; 
		else {
			n_infeasible_gen = 0;
			if(update_best == 1) n_non_updating_gen = 0;
			else n_non_updating_gen++;	
		}

		printf("number of infeasible individuals: %d\n", n_infeasible_ind);
		printf("n_non_updating: %d\n", n_non_updating_gen);
		printf("number of consecutive infeasible generations: %d\n", n_infeasible_gen);

		next_generation(pop1, pop2, fitness_values);

		unsigned char *temp=pop1; pop1=pop2; pop2=temp;

		count++;
	} while (n_infeasible_gen < 10 && n_non_updating_gen < 4 && count < 100);

	if (n_infeasible_gen == 10) return 1;

	return 0;
}
