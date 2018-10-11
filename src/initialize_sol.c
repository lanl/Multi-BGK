#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "initialize_sol.h"
#include "units/unit_data.c"
#include "input.h"

void initialize_sol_inhom(double ***f, int nint, double *int_loc, double *ndens_in, double *v_in, double *T_in, int Nx, double *x, int nspec, int Nv, int order, double **c, double *m, double **n_oned, double ***v_oned, double **T_oned) {

  //accessing multid arrays
  //f - (Nx_rank+2*order) nspec Nv*Nv*Nv
  //int_loc - nint
  //rho/v/T_in - nint nspec
  //x - nv
  //m - nspec
  //c - nspec nv

  int i,j,k,l,s;
  int int_id;
  int vIndex,inputIndex;

  for(l=0;l<Nx;l++) {

    //Figure out the interval
    //If no interval is found it defaults to zero
    int_id = 0; 
    for(j=0;j<nint;j++)
      if(x[l] > int_loc[j])
	int_id = j;
    
    for(s=0;s<nspec;s++) {
      inputIndex = int_id*nspec + s;

      //set moment data
      if(ndens_in[inputIndex] != 0) {
	n_oned[l][s] = ndens_in[inputIndex];
	v_oned[l][s][0] = v_in[inputIndex];
	v_oned[l][s][1] = 0.0;
	v_oned[l][s][2] = 0.0;
	T_oned[l][s] = T_in[inputIndex];	
      }
      else {
	n_oned[l][s] = 0.0;
	v_oned[l][s][0] = 0.0;
	v_oned[l][s][1] = 0.0;
	v_oned[l][s][2] = 0.0;
	T_oned[l][s] = 0.0;
      }

      //set distribution data
      for(i=0;i<Nv;i++) 
	for(j=0;j<Nv;j++)
	  for(k=0;k<Nv;k++) {
	    vIndex = k + Nv*(j + Nv*i);	    

 	    if(ndens_in[inputIndex] != 0.0) {
	      f[l+order][s][vIndex] =  ndens_in[inputIndex]*pow(m[s]/(2.0*M_PI*T_in[inputIndex]/ERG_TO_EV_CGS),1.5)*
		exp(-(0.5*m[s]/(T_in[inputIndex]/ERG_TO_EV_CGS))*
		    ( (c[s][i]-v_in[inputIndex])*(c[s][i]-v_in[inputIndex]) 
	            + (c[s][j]*c[s][j])
		    + (c[s][k]*c[s][k]) ));
	    }
	    else {
	      f[l+order][s][vIndex] = 0.0;
	    }
	  }	
    }
  }


  
}


void initialize_sol_load_inhom_file(int Nx, int nspec, double **n_oned, double ***v_oned, double **T_oned, char *input_file_data_filename) 
{

  /*
    input file layout
  
    Do it by species, e.g. 

    species
    0
    n
    1.0 1.0 1.0 ... 1.0
    v
    1.0 1.0 1.0 ... 1.0
    T
    1.0 1.0 1.0 ... 1.0

    will probably need to use strtok
  */

  // add input prefix

  char input_path[256] = "./input/";

  printf("In load file, ready to load file %s\n", input_file_data_filename);
  
  strcat(input_path,input_file_data_filename);

  printf("Input path: %s\n",input_path);

  FILE *data_file = fopen(input_path,"r");
  
  char line[10000] = {"dummy"};
  char *token;
  int stopFlag = 0;

  int sp; 
  int l;

  while(!stopFlag) {
    read_line(data_file,line);

    printf("%s\n",line);

    if(strcmp(line,"Stop") == 0) 
      stopFlag = 1;
    else if(strcmp(line,"species") == 0) {
      sp = read_int(data_file);

      read_line(data_file,line);
      printf("%s\n",line);      
      if(strcmp(line,"dens") != 0) {
	printf("Error in initial condition data file - expected 'dens'\n");
	exit(1);
      }
      else {
	read_line(data_file,line);
	token = strtok(line," ");	    
	for(l=0;l<Nx;l++) {	  
	  n_oned[l][sp] = atof(token);
	  token = strtok(NULL," ");	    
	}
      }

      read_line(data_file,line);
      printf("%s\n",line);

      if(strcmp(line,"velocity") != 0) {
	printf("Error in initial condition data file - expected 'velocity'\n");
	exit(1);
      }
      else {
	read_line(data_file,line);
	token = strtok(line," ");	    
	for(l=0;l<Nx;l++) {
	  v_oned[l][sp][0] = atof(token);
	  v_oned[l][sp][0] = 0;
	  v_oned[l][sp][0] = 0;
	  token = strtok(NULL," ");	    
	}
      }

      read_line(data_file,line);
      printf("%s\n",line);
      if(strcmp(line,"temperature") != 0) {
	printf("Error in initial condition data file - expected 'temperature'\n");
	exit(1);
      }
      else {
	read_line(data_file,line);
	token = strtok(line," ");	    
	for(l=0;l<Nx;l++) {
	  T_oned[l][sp] = atof(token);
	  token = strtok(NULL," ");	    
	}
      }

    }
    else {
      printf("Error in initial condition data file - missed species block or Stop command\n");
      exit(1);
    }
      
  }

  // check data load
  for(sp=0;sp<nspec;sp++) {
    printf("Species %d \n",sp);
    for(l=0;l<Nx;l++) {
      printf("%d: n: %g v: %g T: %g\n",l,n_oned[l][sp],v_oned[l][sp][0],T_oned[l][sp]);
    }
  }

}
void initialize_sol_inhom_file(double ***f, int Nx, int nspec, int Nv, int order, double **c, double *m, double **n_oned, double ***v_oned, double **T_oned) 
{
  //accessing multid arrays
  //f - (Nx+2*order) nspec Nv*Nv*Nv
  //rho/v/T_in - Nx nspec
  //x - Nx
  //m - nspec
  //c - nspec nv

  int i,j,k,l,s;
  int vIndex;

  for(l=0;l<Nx;l++) {
    for(s=0;s<nspec;s++) {
      //set distribution data
      for(i=0;i<Nv;i++) 
	for(j=0;j<Nv;j++)
	  for(k=0;k<Nv;k++) {
	    vIndex = k + Nv*(j + Nv*i);	    

 	    if(n_oned[l][s] != 0) {
	      f[l+order][s][vIndex] =  n_oned[l][s]*pow(m[s]/(2.0*M_PI*T_oned[l][s]/ERG_TO_EV_CGS),1.5)*
		exp(-(0.5*m[s]/(T_oned[l][s]/ERG_TO_EV_CGS))*
		    ( (c[s][i]-v_oned[l][s][0])*(c[s][i]-v_oned[l][s][0]) 
	            + (c[s][j]*c[s][j])
    	            + (c[s][k]*c[s][k]) ));
	    }
	    else {
	      f[l+order][s][vIndex] = 0.0;
	    }
	  }	
    }
  }
  


}
