#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

typedef struct Network *network;
struct Network {
	unsigned int size;
	struct node *individuals;
};
 
struct node {
	bool state;
	float threshold;
	unsigned int degree;
	unsigned int *connections;
};
 
typedef struct Threshold_distribution threshold_distribution;
struct Threshold_distribution {
	float mean;
	float stdev;
};
 
/* General procedure */
network create_network(int, unsigned int, unsigned int, gsl_rng *);
/* Subprocedures */
void generate_metapopulation(network);
void generate_2d_lattice_rect(network);
void generate_2d_lattice_cross(network);
void generate_regular_random_network(network, gsl_rng *);
/* Termination procedure */
void delete_network(network);

void do_update(network *, gsl_rng *);
void do_swap(network, gsl_rng *); 
void do_swap_in_the_nhood(network, gsl_rng *);

unsigned int total(network);
unsigned int will_go(network);

void shuffle(int *, int, gsl_rng *);
float average(float *, int);
float stdev(float *, int);
float randgauss(float, float, gsl_rng *);
long seedgen(void);
float ran1(long *);