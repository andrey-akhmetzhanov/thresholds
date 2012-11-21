#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
#include <gsl/gsl_rng.h>

/* For getpid() function */
#include <sys/types.h>
#include <unistd.h>

#include "tsp_module.h"

network create_network(int network_type, unsigned int characteristic_dimension, unsigned int characteristic_degree, gsl_rng *rng)
{
	int i,j;
	network nt = (network)malloc(sizeof(struct Network));
	if (network_type==1) 
	{
		nt->size = characteristic_dimension*(characteristic_degree+1);
		nt->individuals = (struct node*)malloc(nt->size*sizeof(struct node)); 
		for (i=0;i<nt->size;i++)
		{
			nt->individuals[i].degree = characteristic_degree-1;
			nt->individuals[i].connections = malloc(nt->individuals[i].degree*sizeof(unsigned int));
		}
		generate_metapopulation(nt);
	}
	else if (network_type==21)
	{		
		nt->size = characteristic_dimension*characteristic_dimension;
		nt->individuals = (struct node*)malloc(nt->size*sizeof(struct node)); 
		for (i=0;i<nt->size;i++)
		{
			nt->individuals[i].degree = (2*characteristic_degree+1)*(2*characteristic_degree+1)-1;
			nt->individuals[i].connections = malloc(nt->individuals[i].degree*sizeof(unsigned int));
		}
		generate_2d_lattice_rect(nt);
	}
	else if (network_type==22)
	{		
		nt->size = characteristic_dimension*characteristic_dimension;
		nt->individuals = (struct node*)malloc(nt->size*sizeof(struct node)); 
		for (i=0;i<nt->size;i++)
		{
			nt->individuals[i].degree = 2*characteristic_degree*(characteristic_degree+1);
			nt->individuals[i].connections = malloc(nt->individuals[i].degree*sizeof(unsigned int));
		}
		generate_2d_lattice_cross(nt);
	}
	else if (network_type==3)
	{		
		nt->size = characteristic_dimension*(characteristic_degree+1);
		nt->individuals = (struct node*)malloc(nt->size*sizeof(struct node)); 
		for (i=0;i<nt->size;i++)
		{
			nt->individuals[i].degree = characteristic_degree;
			nt->individuals[i].connections = malloc(nt->individuals[i].degree*sizeof(unsigned int));
		}
		generate_regular_random_network(nt,rng);
	}
	else 
	{
		fprintf(stderr,"Error: unrecognized type of the network.\n");
		exit(1);
	}

	return nt;
}

void generate_metapopulation(network nt)
{
	int size_subpopulation = nt->individuals[0].degree+1;
	for (int i=0;i<nt->size;i++)
	{
		int k = 0;
		for (int j=0;j<size_subpopulation;j++)
			if (j!=i%size_subpopulation) 
				nt->individuals[i].connections[k++] = j+size_subpopulation*(i/size_subpopulation);
	}
}

void generate_2d_lattice_rect(network nt)
{
	int i,j,k,ii,jj;
	int R = (sqrt(nt->individuals[0].degree+1)-1)/2, L = sqrt(nt->size);
	for (i=0;i<nt->size;i++)
	{
		k = 0;
		for (ii=-R;ii<=R;ii++)
			for (jj=-R;jj<=R;jj++)
				if (!((ii==0)&&(jj==0)))
					nt->individuals[i].connections[k++] = ((i%L+ii)%L+L)%L+L*(((i/L+jj)%L+L)%L);
	}
}

void generate_2d_lattice_cross(network nt)
{
	int i,j,k,ii,jj;
	int R = (sqrt(2*nt->individuals[0].degree+1)-1)/2, L = sqrt(nt->size);
	for (i=0;i<nt->size;i++)
	{
		k = 0;
		for (ii=-R;ii<=R;ii++)
			for (jj=-R;jj<=R;jj++)
				if ((!((ii==0)&&(jj==0)))&&(ii*ii+jj*jj<=R*R))
					nt->individuals[i].connections[k++] = ((i%L+ii)%L+L)%L+L*(((i/L+jj)%L+L)%L);
	}
}

void generate_regular_random_network(network nt, gsl_rng *rng)
{
	/* Generation of the network according to Newman */
	int i,j,k;
	unsigned int total_number_of_connections = 0;
	for (i=0;i<nt->size;i++)
		total_number_of_connections += nt->individuals[i].degree;
	int *set, last_element = 0; 
	set = malloc(total_number_of_connections*sizeof(unsigned int));
	int last = 0;
	for (i=0;i<nt->size;i++)
		for (j = 0;j<nt->individuals[i].degree;j++)
		{
			*(set+last) = i; last++;
		} 
	shuffle(set,total_number_of_connections,rng);
 	int clast[nt->size];
	for (i=0;i<nt->size;i++) clast[i] = 0;
	for (i=0;i<(total_number_of_connections-total_number_of_connections%2)/2;i++) 
	{
		if (set[2*i]!=set[2*i+1])
		{
			nt->individuals[set[2*i]].connections[clast[set[2*i]]++] = set[2*i+1];
			nt->individuals[set[2*i+1]].connections[clast[set[2*i+1]]++] = set[2*i];
		}
		else 
			nt->individuals[set[2*i]].degree--;
	}
	if (total_number_of_connections%2!=0)
		nt->individuals[set[total_number_of_connections]].degree--;
	free(set);	
}

void delete_network(network nt)
{
	for (int i=0; i<nt->size; i++)
		free(nt->individuals[i].connections);
	free(nt->individuals);
	free(nt);
}

/* Update the state of randomly picked individual */
void do_update(network *nt, gsl_rng *rng)
{
	int i=gsl_rng_uniform(rng)*(*nt)->size, j, sum_local=0;
	for (j=0;j<(*nt)->individuals[i].degree;j++) 
		sum_local+=(*nt)->individuals[(*nt)->individuals[i].connections[j]].state;
	if (sum_local>=((*nt)->individuals[i].threshold)*((*nt)->individuals[i].degree))
		(*nt)->individuals[i].state = 1;
	else
		(*nt)->individuals[i].state = 0;
}

/* Swap two randomly picked individuals */
void do_swap(network nt, gsl_rng *rng)
{
	int i = gsl_rng_uniform(rng)*nt->size, j = gsl_rng_uniform(rng)*nt->size; 
	bool st_tmp, *st1, *st2;
	st1 = &(nt->individuals[i].state); st2 = &(nt->individuals[j].state);
	st_tmp = *st1; *st1 = *st2; *st2 = st_tmp;
}

void do_swap_in_the_nhood(network nt,gsl_rng *rng)
{
	int i = gsl_rng_uniform(rng)*nt->size, j = nt->individuals[i].connections[(int)(gsl_rng_uniform(rng)*nt->individuals[i].degree)];
	bool st_tmp, *st1, *st2;
	st1 = &(nt->individuals[i].state); st2 = &(nt->individuals[j].state);
	st_tmp = *st1; *st1 = *st2; *st2 = st_tmp;
}

unsigned int total(network nt)
{
	unsigned int i, sum = 0;
	for (i=0;i<nt->size;i++)
		sum += nt->individuals[i].state;
	return(sum);
}
 
unsigned int will_go(network nt)
{
	unsigned int i, j, sum_local, sum = 0;
	for (i=0;i<nt->size;i++)
	{
		sum_local = 0; j = 0;
		for (j=0;j<nt->individuals[i].degree;j++) 
			sum_local += nt->individuals[nt->individuals[i].connections[j]].state;
		while (j<nt->individuals[i].degree);
		if (sum_local>=nt->individuals[i].threshold*nt->individuals[i].degree)
			sum += 1;
	}
	return(sum);
}
 
void shuffle(int *array, int n, gsl_rng *rng) 
{
  	int i, j, tmp;
  	for (i=n-1;i>0;i--) 
  	{
    	j = (int)(gsl_rng_uniform(rng)*(i+1));
    	tmp = array[j];
    	array[j] = array[i];
    	array[i] = tmp;
  	}
}

float average(float *array, int n) 
{
	int i;
	float sum = 0.0;
	for (i=0;i<n;i++)
    	sum += array[i];
    return(sum/n);
}

float stdev(float *array, int n) 
{
	float aver = average(array,n);
	int i;
	float sum = 0.0;
	for (i=0;i<n;i++) 
    	sum += 1.0*(array[i]-aver)*(array[i]-aver);
    return(sqrt(sum/(n-1)));
}

float randgauss(float mu, float sigma, gsl_rng *rng)
{
	float x1, x2, w;
	do {
		x1 = 2.0*gsl_rng_uniform(rng)-1.0;
		x2 = 2.0*gsl_rng_uniform(rng)-1.0;
		w = x1*x1 + x2*x2;
	} while (w >= 1.0);
	w = sqrt((-2.0*log(w))/w);
	return(sigma*x1*w+mu);
}

long seedgen(void)  {
/* From section 7.1 in Random Numbers In Scientific Computing: An Introduction by Katzgrabber
 * This algorithm is suggested to use for parallel programming. It gives a hash function of time
 * and PID to generate a seed. */

    long s, seed, pid;
    time_t seconds;

    pid = getpid();
    s = time(&seconds); /* get CPU seconds since 01/01/1970 */

    seed = abs(((s*181)*((pid-83)*359))%104729); 
    return seed;
}