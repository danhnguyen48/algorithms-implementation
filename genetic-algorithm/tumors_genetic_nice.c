#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "RKF78.h"
#include "rand1.h"

#define lambda_par 			 0.336
#define lambda_log_theta_par 9.284023094951992774680543277273261947643
#define NZero_par 			 20000
#define pop_size 			 30000
#define n_par 				 11		//time steps
#define d_par 	 			 10		//n of treatments in a cure
#define m_par 				 4		//organs that cannot be damaged
#define CMAXj_par	 		 15		
#define CCUMj_par	 		 127	//constraint 2
#define NMax_par 			 9.5e+11//max number of tumor cells
#define MAXDOUBLE			 1.9e+16
#define prob_cross			 0.3	//probability of crossover occurrence
#define prob_mut			 0.1	//proability of mutation occurrence
#define tournament_size      10		//number of selections with replacement in each choice of parents

static float k_j[d_par] = {0.12, 0.0502, 0.0637, 0.1347, 0.0902, 0.0546, 0.0767, 0.1121, 0.0971, 0.0403};
static float t_i[n_par] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 33.0};
static float eta_kj[m_par][d_par] = { //curative parameters 
	{0.0036, 0.0098, 0.0061, 0.0009, 0.0003, 0.0108, 0.0045, 0.0021, 0.0096, 0.0125},
	{0.0063, 0.0082, 0.0062, 0.0062, 0.0083, 0.013, 0.0039, 0.0019, 0.0015, 0.005},
	{0.0129, 0.0018, 0.0116, 0.0021, 0.009, 0.0129, 0.0054, 0.0049, 0.0093, 0.0066},
	{0.0053, 0.0086, 0.0067, 0.0029, 0.0089, 0.0054, 0.0042, 0.0095, 0.0112, 0.0092} }; 
/* //palliative parameters
{
	{0.00612, 0.01666, 0.01037, 0.00153, 0.00051, 0.01836, 0.00765, 0.00357, 0.01632, 0.02125},
	{0.01071, 0.01394, 0.01054, 0.01054, 0.01411, 0.0221, 0.00663, 0.00323, 0.00255, 0.0085},
	{0.02193, 0.00306, 0.01972, 0.00357, 0.0153, 0.02193, 0.00918, 0.00833, 0.01581, 0.01122},
	{0.00901, 0.01462, 0.01139, 0.00493, 0.01513, 0.00918, 0.00714, 0.01615, 0.01904, 0.01564}}; */
	
typedef struct { double drift_i; } ODE_Parameters;


void gompertz(double t, double N, double *der, void *Params){

	*der = ((N < 1.e-16) ? 0.0 : N*(lambda_log_theta_par - lambda_par*log(N) - ((ODE_Parameters *) Params)->drift_i));
}

void mod_gompertz(double t, double N, double *der, void* Params){
	*der = ((N < 1.e-16) ? 0.0 : N*(lambda_log_theta_par - lambda_par*log(N)));
}

int init_population(unsigned char* pop){
	//fills the memory pointed by pop with individuals 
	unsigned char* drive=pop;
	long int seed=-rand();
	int generated=0;
	while (drive != pop+pop_size*d_par*(n_par-1)){

		//generating random individual
		for (int i=0; i<d_par*(n_par-1); i++){
			generated=(int) (ran1(&seed)*1.0e+4);
			*(drive+i)=generated % CMAXj_par;
		}

		drive+=d_par*(n_par-1);
	}

	return 0;
}

unsigned char test_if_constraints_2and_3_are_verified(unsigned char *Cij){
	register unsigned char i, j, k;
	unsigned char npar = n_par - 1, *Cijofi;

	//constraint 2
	for(j=0; j < d_par; j++){ 
		unsigned int ccumj = 0U;
		for(i=0; i < npar; i++) ccumj += *(Cij + i*d_par + j);
		if(ccumj > CCUMj_par) return 0U;
	}

	//constraint 3
	for(i=0, Cijofi=Cij; i < npar; i++, Cijofi += d_par) {
		for(k=0; k < m_par; k++){ 
			double Cseffk = 0.0;
			for(j=0; j < d_par; j++) Cseffk += eta_kj[k][j] * Cijofi[j];
			if(Cseffk > 1.0) return 0U;
		}
	}
	return 1U;
}

double curative_fitness(unsigned char *Cij){
	register unsigned char i, j;
	ODE_Parameters GompertzParams;
	double N = NZero_par, t = t_i[0];
	double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
	unsigned char curativecounter = 0U, npar = n_par - 1;
	double integral = 0.0, lastt = t, lastN = N;

	if(!test_if_constraints_2and_3_are_verified(Cij)) return MAXDOUBLE;

	

	for(i=0; i < npar; i++){ 
		double tfin = t_i[i+1]; // implementing treatment i
		GompertzParams.drift_i = 0.0;
		for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

		while(t+h < tfin) {
			RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
			if(N > NMax_par) return MAXDOUBLE;
			integral += (lastN + N)*(t - lastt);
			lastt = t; lastN = N;
		}

		do { 
		h = tfin - t;
		RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
		if(N > NMax_par) return MAXDOUBLE;
		integral += (lastN + N)*(t - lastt);
		lastt = t; lastN = N;
		} while (t < tfin);

		if(N < 1000) { 
			curativecounter++; 
			if(curativecounter > 2) return integral/2.0; 
		} else curativecounter = 0U;
	}

	return MAXDOUBLE;
}

double palliative_fitness(unsigned char* Cij){
	/*returns a measure of 1/lifespan, function that we want to minimize. If the patient dies before the end of 
	  the treatment we consider that 1/lifespan is infinity*/
	register unsigned char i, j;
	ODE_Parameters GompertzParams;
	double N = NZero_par, t = t_i[0];
	double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
	unsigned char npar = n_par - 1;
	double integral = 0.0, lastt = t; 

	if(!test_if_constraints_2and_3_are_verified(Cij)) return MAXDOUBLE;

	//see if the patient survives until the end of the treatment
	for(i=0; i < npar; i++){ 
		double tfin = t_i[i+1]; // Implementing treatment i
		GompertzParams.drift_i = 0.0;
		for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

		while(t+h < tfin) {
			RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
			if(N > NMax_par) return MAXDOUBLE;
			lastt = t;
		}

		do { 
		h = tfin - t;
		RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, gompertz);
		if(N > NMax_par) return MAXDOUBLE;
		lastt = t; 
		} while (t < tfin);

	}

	//if the patient survives calculate the effective lifespan
	while(N<NMax_par){
		RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, mod_gompertz);
		lastt=t;
	}

	return 1/lastt;
}

int select_with_replacement (char* pop, double* fitness, unsigned char* parent){
	//selects an individual from the population using the tournament selection principle
	
	long int seed=-rand();
	int selection=((int) (ran1(&seed)*1.0e+4))%pop_size;
	unsigned char* best=pop+selection*d_par*(n_par-1);
	double fit_best=fitness[selection];
	
	unsigned char* next=NULL;
	for (int i=2; i<=tournament_size; i++){
		//new individual
		selection=((int)(ran1(&seed)*1.0e+4))%pop_size;
		next= pop+selection*d_par*(n_par-1);
		//if its fitness is better than the best one it becomes the best
		if (fitness[selection]<fit_best){
			fit_best=fitness[selection]; 
			best=next;}
	}

	//copying the selected individual
	for (unsigned char i=0; i<d_par*(n_par-1); i++){
		*(parent+i)=*(best+i);
	}
	
	return 1;
}

int crossover(unsigned char* parent_a, unsigned char* parent_b, unsigned char* child1, unsigned char* child2){
	//uniform crossover between the parents generates the children
	long int seed=-rand();
	double p=0;

	for (int i=0; i<d_par*(n_par-1); i++){
		p=ran1(&seed);
		
		if (p<=prob_cross){ //no swapping 
			child1[i]=parent_a[i];
			child2[i]=parent_b[i];
			
		} else { //swapping
			child1[i]=parent_b[i];
			child2[i]=parent_a[i];
			
		}
	}
	
	return 1;

}

int mutate(unsigned char* child){
	//bitflip mutation
	long int seed=-rand();
	double p=0;
	for (int i=0; i<d_par*(n_par-1); i++){
		p=ran1(&seed);
		if (p<= prob_mut) child[i]=CMAXj_par-child[i];
	}
	return 1;
}

int next_generation(unsigned char* pop1, unsigned char* pop2, double* fitness_value){
			unsigned char* parent_a=malloc(d_par*(n_par-1)*sizeof(unsigned char));
			unsigned char* parent_b=malloc(d_par*(n_par-1)*sizeof(unsigned char));
			unsigned char* child1=malloc(d_par*(n_par-1)*sizeof(unsigned char));
			unsigned char* child2=malloc(d_par*(n_par-1)*sizeof(unsigned char));
			for (int i=0; i<pop_size/2; i++){
				
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
					pop2[i*d_par*(n_par-1)+j]=child1[j];
					pop2[pop_size/2+i*d_par*(n_par-1)+j]=child2[j];
				}	
			}
			//free memory 
			free(parent_a); free(parent_b); free(child1); free(child2);
}



int look_for_curative(unsigned char* pop1, unsigned char* pop2, double* best_fitness, unsigned char* best, double* fitness_values){
	//LOOKING FOR A CURATIVE SOLUTION
	//control variables
	int n_infis_gen=0;			//number of consecutive infeasible generations
	int n_infis_ind=0;			//number of infeasible individuals in a generation
	int n_non_updating_gen=0; 	//number of consecutive generations in which the best individual does not change
	int update_best;			//belongs to {0,1}: indicates if the best indiv. has been updated in the current generation
	int best_pos=0;				//position of the best element of the pop in the vector

	int gen_counter=0;
	do {
		printf("GENERATION %d\n", gen_counter);
		n_infis_ind=0; //reset the number of infeasibles
		update_best=0; 
		
		//initializing the fitness values, computing the best position, counting the infeasibles
		for(int i=0; i<pop_size; i++){
			fitness_values[i]=curative_fitness(pop1+i*d_par*(n_par-1));
			
			if (fitness_values[i]==MAXDOUBLE) n_infis_ind++; 
			else if (fitness_values[i]<*best_fitness){ 

					best_pos=i;
					*best_fitness=fitness_values[i];
					update_best=1;
				}
			}

		//saving the best individual if it has changed in the generation
		if (update_best==1){for (unsigned char j=0; j<d_par*(n_par-1); j++){ *(best+j)=*(pop1+best_pos*d_par*(n_par-1)+j);}}
			
		//if the generation is infeasible increase n_infis_gen, otherwise set it to 0 and check if the best is updated
		if (n_infis_ind==pop_size) n_infis_gen++; else {
			n_infis_gen=0;
			//if the best is updated and the generation is not infeasible the number of non-updating gen is set to 0, otherwise it increases
			if(update_best==1) n_non_updating_gen=0; else n_non_updating_gen++;
		}

		printf("%d infeasible individuals, ", n_infis_ind);
		printf("n_non_updating= %d, ", n_non_updating_gen);
		printf("n_inf_generations= %d\n", n_infis_gen);

		//creating the next generation
		next_generation(pop1, pop2, fitness_values);

		//switching the pointers of pop1 and pop2
		unsigned char* temp=pop1; pop1=pop2; pop2=temp;

		gen_counter++;
	} while (n_infis_gen<10 & n_non_updating_gen<4 & gen_counter<100);

	if (n_infis_gen==10) return 1;

	return 0;
}

int look_for_palliative(unsigned char* pop1, unsigned char* pop2, double* best_fitness, unsigned char* best, double* fitness_values){
	//LOOK FOR A PALLIATIVE SOLUTION
	//control variables
	int n_infis_gen=0;			//number of consecutive infeasible generations
	int n_infis_ind=0;			//number of infeasible individuals in a generation
	int n_non_updating_gen=0; 	//number of consecutive generations in which the best individual does not change
	int update_best;			//belongs to {0,1}: indicates if the best indiv. has been updated in the current generation
	int best_pos=0;				//position of the best element of the pop in the vector

	int gen_counter=0;
	do {
		printf("GENERATION %d\n", gen_counter);
		n_infis_ind=0;
		update_best=0;

		//initializing the fitness values, computing the best position, counting the infeasibles
		for (int i=0; i<pop_size; i++){
			fitness_values[i]=palliative_fitness(pop1+i*d_par*(n_par-1)); 

			if (fitness_values[i]==MAXDOUBLE) {n_infis_ind++;} 
			else if (fitness_values[i]<*best_fitness){

				best_pos=i;
				*best_fitness=fitness_values[i];
				update_best=1;
			}
		}

		//saving the best individual if it has changed in the generation
		if(update_best==1) {for (int j=0; j<d_par*(n_par-1); j++){ *(best+j)=*(pop1+best_pos*d_par*(n_par-1)+j);}}
			
		//if the generation is infeasible increase n_infis_gen, otherwise set it to 0 and check if the best is updated
		if (n_infis_ind==pop_size) n_infis_gen++; else {
			n_infis_gen=0;
			//if the best is updated and the generation is not infeasible non-updating gen is set to 0, otherwise it increases
			if(update_best==1) n_non_updating_gen=0; else n_non_updating_gen++;	
		}

		printf("%d infeasible individuals, ", n_infis_ind);
		printf("n_non_updating= %d, ", n_non_updating_gen);
		printf("n_inf_generations= %d\n", n_infis_gen);
				
		//creating the next generation
		next_generation(pop1, pop2, fitness_values);

		//switching the pointers of pop1 and pop2
		unsigned char* temp=pop1; pop1=pop2; pop2=temp;

		gen_counter++;
	} while (n_infis_gen<10 & n_non_updating_gen<4 & gen_counter<100);

	if (n_infis_gen==10) return 1;

	return 0;
}

int main(int argc, char *argv[]){
	//memory reservation
	unsigned char* pop1=(unsigned char*) malloc(pop_size*d_par*(n_par-1)*sizeof(unsigned char));
	unsigned char* pop2=(unsigned char*) malloc(pop_size*d_par*(n_par-1)*sizeof(unsigned char));
	double* fitness_values=(double*) malloc(pop_size*sizeof(double));			//contains the ordered fitness values of a generation
	
	unsigned char* best=(unsigned char*) malloc(d_par*(n_par-1)*sizeof(unsigned char));	//best individual founded by the algorithm
	double best_fitness=MAXDOUBLE;
	
	init_population(pop1);	//initialization of the population

	//look for a curative solution
	int curative=look_for_curative(pop1, pop2, &best_fitness, best, fitness_values);

	if (curative==0) {

		printf("A curative treatment has been found. Its doses are \n");
		for (int i=1; i<=d_par*(n_par-1); i++) {
			printf("%d ",best[i-1]);
			if (i%d_par==0) printf("\n");
		}
		printf("and its fitness value is %f ", best_fitness);
		printf("\n");

	} else { 
		//look for a palliative solution
		best_fitness=MAXDOUBLE; //reset the best fitness value

		int palliative=look_for_palliative(pop1, pop2, &best_fitness, best, fitness_values);

		if(palliative==0){

			printf("A palliative treatment has been found. Its doses are \n");
			for (int i=1; i<=d_par*(n_par-1); i++) {
				printf("%d ",best[i-1]);
				if (i%d_par==0) printf("\n");
			}
			printf("and the lifespan of the patient is %f ", 1/best_fitness);
			printf("\n");

		} else printf("Neither a curative nor a palliative treatment have been found");
	}
	

}