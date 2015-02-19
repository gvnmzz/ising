/*    University of Warwick -- CO904 Statistical Mechanics 2013    */
/* Ising model simulation based on code by Charo del Genio -- 2012 */

/* compile using: gcc Ising2d.c -o Ising2d */
/* run using: ./Ising2d [size] [initprob] [temperature] [number of sweeps] [seed] [output filename] */

/* BEGIN HEADER */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* global variables */
unsigned int STATE[32];
unsigned long int seed;

/* prototypes */
void parameters(int nump, char *param[], int *size, double *initprob, double *temperature, int *sweeps, FILE **p_f_output);
inline int  DeltaIndex(const int x, const int y, const int size, int **spin);
/* END HEADER */

//=============================
int main(int argc, char *argv[])
{
  char outputfile[80];
  int size, sweeps, DeltaE, x, y, diff, square, **spin;
  register int i, j, E=0, M=0;
  double initprob, temperature, Factors[5];
  FILE *p_f_output;
  FILE *spin_output;

  // computing |m| during simulation (bad practice -- demonstration only)
  double abs_av_m=0.0;  
  int t_av;

  /* Check parameters */
  parameters(argc, argv, &size, &initprob, &temperature, &sweeps, &p_f_output);
  
  square = size*size; // square lattice only
  spin = malloc(size*sizeof(int*)); // allocate lattice
  spin[0] = malloc(square*sizeof(int));
  for (i=1; i<size; i++) spin[i] = spin[i-1] + size;
  
  for (i=0; i<size; i++) for (j=0; j<size; j++) spin[i][j] = initprob > drand48() ? 1 : -1; // Create initial state
  
  /* Boltzmann factors */
  Factors[0]=exp(-8.0/temperature); // All spins are equal to the one seeded	DeltaE =  8
  Factors[1]=exp(-4.0/temperature); // 1 spin  is  different			DeltaE =  4
  Factors[2]=1.0;		    // 2 spins are different			DeltaE =  0
  Factors[3]=exp(4.0/temperature);  // 3 spins are different			DeltaE = -4
  Factors[4]=exp(8.0/temperature);  // All spins are different			DeltaE = -8
  
  /* record initial state */
  for (x=0; x<size; x++) for (y=0; y<size; y++) M += spin[x][y];	// Initial value of M
  
  /* print to file |m| */ // initial state
  //  fprintf(p_f_output, "%.10f\t",fabs(((double) M)/square));
  
  for (x=0; x<size; x++) for (y=0; y<size; y++) {	// Initial value of E
      E -= spin[x][y] * spin[x==size-1 ? 0 : x+1][y];
      E -= spin[x][y] * spin[x][y==size-1 ? 0 : y+1];
    }

  /* print to file E */ // initial state
  //  fprintf(p_f_output, "%d\n",E);
  
  /* Simulate using Metropolis */
  for (i=1; i<=sweeps; i++) {  // for each sweep
    
    for (j=1; j<=square; j++) { // run over the lattice size
      x = size * drand48();     // choose a spin at random
      y = size * drand48();
      
      diff = DeltaIndex(x,y,size,spin);   // look at the neighbours
      if ( drand48() < Factors[diff] ) {  // if it flips
	spin[x][y] *= -1;                 // update it
	M += (2*spin[x][y]);	          // update M
	E += 4*(2-diff);	          // update E
      }
      
    } // end run over lattice size 
    
    /* print to file |m| */
    // calculate the average magnetisation (absolute value) -- bad practice, demonstration only
    if (i>10000) { // discard initialisation -- equilibration time
      t_av = i-10000;
      abs_av_m = abs_av_m * (((double) (t_av-1))/((double) t_av)) + fabs(((double) M)/square)/((double) t_av);
    }
  } // end sweep

  fprintf(p_f_output, "%.10f\t%.10f\n",temperature,abs_av_m);
  printf("Average |m| is: %.10f\n",abs_av_m);

  // output the final spin configuration for visualisation
  spin_output = fopen("spins.dat","w");
  for (x=0; x<size; x++) { 
    for (y=0; y<size; y++) { 
      fprintf(spin_output, "%d\t",(spin[x][y] + 1)/2);
    }
    fprintf(spin_output,"\n"); 
  }
  fclose(spin_output); // close file  
  
  free(spin[0]); // free memory
  free(spin);
  
  fclose(p_f_output); // close file
  
  return EXIT_SUCCESS;
} // end main
//==============================


/* Count how many spins point in a direction different from the one the spin at x,y points towards. Inlined for efficiency */
inline int DeltaIndex(const int x, const int y, const int size, int **spin)
{
  int centre, diff=0;
  
  centre = spin[x][y];
  if ( spin[x==0 ? size-1 : x-1][y] != centre ) diff++;
  if ( spin[x==size-1 ? 0 : x+1][y] != centre ) diff++;
  if ( spin[x][y==0 ? size-1 : y-1] != centre ) diff++;
  if ( spin[x][y==size-1 ? 0 : y+1] != centre ) diff++;
  
  return diff;
} // end DeltaIndex


/* Check user-provided parameters */
void parameters(int nump, char *param[], int *size, double *initprob, double *temperature, int *sweeps, FILE **p_f_output)
{	
  if (nump<7) {
    printf("\nUsage: %s [size] [prob] [temp] [sweeps] [seed] [output]\n\n",param[0]);
    printf("\t[size]   is the dimension of each side of the 2d Ising lattice.\n");
    printf("\t[prob]   is the probability for each spin to initially point up.\n");
    printf("\t[temp]   is the temperature to simulate.\n");
    printf("\t[sweeps] is the number of sweeps to run.\n");
    printf("\t[seed]   is the random number generator seed (positive integer).\n");
    printf("\t[output] is the output file name.\n\n");
    exit(EXIT_FAILURE);
  }
  
  *size = atoi(param[1]);
  if (*size<2) {
    printf("\nThe size of the model must be at least 2x2.\n\n");
    exit(EXIT_FAILURE);
  }
  
  *initprob = atof(param[2]);
  if (*initprob<0.0 || *initprob>1.0) {
    printf("\nThe probability for each spin to initially point up must be between 0 and 1.\n\n");
    exit(EXIT_FAILURE);
  }
  
  *temperature = atof(param[3]);
  
  *sweeps = atoi(param[4]);
  if (*sweeps<1) {
    printf("\nYou must run at least 1 sweep.\n\n");
    exit(EXIT_FAILURE);
  }
  
  seed=strtoul(param[5],NULL,10);
  if (seed<=0) {
    printf("\nThe random number seed must be a positive integer.\n\n");
    exit(EXIT_FAILURE);
  }
  //  seedrng();
  srand48(seed);

  /*    
  if ( (*p_f_output = fopen(param[6],"r") ) != NULL ) {
    printf("\nError!\nOutput file %s already exists!\n\n",param[6]);
    exit(EXIT_FAILURE);
  }
  if ( (*p_f_output = fopen(param[6],"w") ) == NULL ) {
    printf("\nError creating output file %s.\n\n",param[6]);
    exit(EXIT_FAILURE);
  }
  */
  if ( (*p_f_output = fopen(param[6],"a") ) == NULL ) { // currently I want to append to the file
    printf("\nError creating output file %s.\n\n",param[6]);
    exit(EXIT_FAILURE);
  }

  return;
} // end parameters
