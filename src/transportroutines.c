#include <math.h>
#include "transportroutines.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

//Computes the transport term 

static int N;
static int nX;
static int ns;
static int order;
static double Lx;
static double L_v;
static double h_v;
static double dt;
static int ICChoice;
static double TWall;
static double VWall;
static double T0,T1;
static double V0,V1;
static double *x, *dx, **c;
static double ***f_star;

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
double min(double in1, double in2)
{
  return in1 < in2 ? in1 : in2;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
double max(double in1, double in2)
{
  return in1 > in2 ? in1 : in2;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
double minmod(double in1, double in2, double in3)
{  
  if ( (in1 > 0) && (in2 > 0) && (in3  > 0)) {
    return min(min(in1, in2), in3);
  }
  else if ( (in1 < 0) && (in2 < 0) && (in3 < 0)) {
    return max(max(in1, in2), in3);
  }
  else
    return 0;
}

void initialize_transport(int numV, int numX, int nspec, double *xnodes, double *dxnodes, double Lx_val, double **vel, int ord, double timestep) {
  N = numV;
  nX = numX; //Remember that this is nX in the rank
  x = xnodes;
  dx = dxnodes;
  ns = nspec;
  Lx = Lx_val; 

  c = vel;

  order = ord;
  dt = timestep;

  int i,l;
  f_star = malloc((nX+2*order)*sizeof(double *));
  for(l=0;l<nX+2*order;l++){
    f_star[l] = malloc(nspec*sizeof(double *));
    for(i=0;i<nspec;i++)
      f_star[l][i] = malloc(N*N*N*sizeof(double));
  }
  
  x = xnodes;
  dx = dxnodes;
}


/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
void fillGhostCellsPeriodic_firstorder(double ***f, int sp) {

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numRanks);
  MPI_Status status;


  if(numRanks != 1) { 
    printf("ERROR: code does not allow for multiple processes at present!\n");
    exit(1);
  }

  //FILL INTERIOR GHOST CELLS

  if((rank % 2) == 0) { //Have even nodes send first
    
    if(rank != (numRanks - 1)) // Send right data to odd nodes
      MPI_Send(f[nX][sp],N*N*N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);    

    if(rank != 0) // Receive left data from odd nodes
      MPI_Recv(f[0][sp],N*N*N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status);


    if(rank != 0) //Send left data to odd nodes
      MPI_Send(f[1][sp],N*N*N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);

    if(rank != (numRanks - 1)) //Recieve right data from odd nodes
      MPI_Recv(f[nX+1][sp],N*N*N,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);

  }
  else { // Odd nodes recieve first
    MPI_Recv(f[0][sp],N*N*N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status); //Receive left data from even nodes

    if(rank != (numRanks-1))
      MPI_Send(f[nX][sp],N*N*N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD); //Send right data to even nodes

  
    if(rank != (numRanks-1)) // Get right data from even nodes
      MPI_Recv(f[nX+1][sp],N*N*N,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);
    
    if(rank != 0) // Send left data to evens
      MPI_Send(f[1][sp],N*N*N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD); 
  }


  //Now deal with boundary 
  if(rank == (numRanks-1)) { //Rightmost rank sends to then receives from rank 0
    MPI_Send(f[nX][sp],N*N*N,MPI_DOUBLE,0,2,MPI_COMM_WORLD);

    MPI_Recv(f[nX+1][sp],N*N*N,MPI_DOUBLE,numRanks-1,2,MPI_COMM_WORLD,&status);
  }

  if(rank == 0) { //Leftmost rank receives then sends to rank N-1
    MPI_Recv(f[0][sp],N*N*N,MPI_DOUBLE,numRanks-1,2,MPI_COMM_WORLD,&status);

    MPI_Send(f[1][sp],N*N*N,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
  }
}

void fillGhostCellsPeriodic_secondorder(double ***f, int sp) {

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numRanks);
  MPI_Status status;


  if(numRanks != 1) { 
    printf("ERROR: code does not allow for multiple processes at present!\n");
    exit(1);
  }

  //FILL INTERIOR GHOST CELLS

  if((rank % 2) == 0) { //Have even nodes send first
    
    if(rank != (numRanks - 1)) {// Send right data to odd nodes
      MPI_Send(f[nX][sp],N*N*N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);    
      MPI_Send(f[nX+1][sp],N*N*N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);    
    }
    if(rank != 0) {// Receive left data from odd nodes
      MPI_Recv(f[1][sp],N*N*N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status);
      MPI_Recv(f[0][sp],N*N*N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status);
    }

    if(rank != 0) {//Send left data to odd nodes
      MPI_Send(f[2][sp],N*N*N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
      MPI_Send(f[3][sp],N*N*N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
    }
    if(rank != (numRanks - 1)) {//Recieve right data from odd nodes
      MPI_Recv(f[nX+3][sp],N*N*N,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);
      MPI_Recv(f[nX+2][sp],N*N*N,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);
    }
  }
  else { // Odd nodes recieve first
    MPI_Recv(f[1][sp],N*N*N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status); //Receive left data from even nodes
    MPI_Recv(f[0][sp],N*N*N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status); //Receive left data from even nodes

    if(rank != (numRanks-1)) {
      MPI_Send(f[nX][sp],N*N*N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD); //Send right data to even nodes
      MPI_Send(f[nX+1][sp],N*N*N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD); //Send right data to even nodes
    }
  
    if(rank != (numRanks-1)) {// Get right data from even nodes
      MPI_Recv(f[nX+3][sp],N*N*N,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);
      MPI_Recv(f[nX+2][sp],N*N*N,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);
    }
    if(rank != 0) {// Send left data to evens
      MPI_Send(f[2][sp],N*N*N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD); 
      MPI_Send(f[3][sp],N*N*N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD); 
    }
  }


  //Now deal with boundary 
  if(rank == (numRanks-1)) { //Rightmost rank sends to then receives from rank 0
    MPI_Send(f[nX],N*N*N,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    MPI_Send(f[nX+1],N*N*N,MPI_DOUBLE,0,2,MPI_COMM_WORLD);

    MPI_Recv(f[nX+3],N*N*N,MPI_DOUBLE,numRanks-1,2,MPI_COMM_WORLD,&status);
    MPI_Recv(f[nX+2],N*N*N,MPI_DOUBLE,numRanks-1,2,MPI_COMM_WORLD,&status);
  }

  if(rank == 0) { //Leftmost rank receives then sends to rank N-1
    MPI_Recv(f[1],N*N*N,MPI_DOUBLE,numRanks-1,2,MPI_COMM_WORLD,&status);
    MPI_Recv(f[0],N*N*N,MPI_DOUBLE,numRanks-1,2,MPI_COMM_WORLD,&status);

    MPI_Send(f[2],N*N*N,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    MPI_Send(f[3],N*N*N,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
  }
}




//Computes the x direction first order upwind solution FOR A SINGLE SPECIES
void upwindOne_x(double ***f, double ***f_conv, double *v, int sp) {
  int i,j,k,l;
  int index;
  double CFL_NUM, CFL_NUM2;

  fillGhostCellsPeriodic_firstorder(f,sp);

  //Note - the 'real' data lives in 1...Nx, the ghost points are 0 and Nx+1

  for(l=1;l<nX+1;l++)  
    //#pragma omp parallel for private(CFL_NUM,i,j,k,index)
    for(i=0;i<N;i++) {
      CFL_NUM = dt*v[i]/dx[l];
      for(j=0;j<N;j++) 
	for(k=0;k<N;k++) {	
	  index = k + N*(j + N*i);
	  //the upwinding
	  if(i < N/2) 
	    f_conv[l][sp][index] = f[l][sp][index] - CFL_NUM*(f[l+1][sp][index] - f[l][sp][index]);	  
	  else 
	    f_conv[l][sp][index] = f[l][sp][index] - CFL_NUM*(f[l][sp][index]   - f[l-1][sp][index]);
	}
    }
      
}

//Computes the v direction first order upwind solution FOR A SINGLE SPECIES
void upwindOne_v(double ***f, double ***f_conv, double *PoisPot, double **qm, double m, double *v, int sp) {
  int i,j,k,l;
  int index;
  double CFL_NUM;
  h_v = v[1] - v[0];

  int rank, numNodes;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numNodes);
  MPI_Status status;
  //int numamt;

  //Set electric field acceleration term
  //This is given by
  // - Z (e_0 E) / m
  //
  // e_0 E = (e phi)_x is given in units of ergs
  // Z is the partial ionization, provided as qm in the input (no unit)
  // m is the mass of the species in grams
  //
  //Thus the acceleration term units are cm / s^2

  double *E = malloc(nX*sizeof(double));
  E[0] = -qm[0][sp]/m*(PoisPot[1]-PoisPot[nX-1])/(2.0*dx[0]);
  for(i=1;i<nX-1;i++) {
    E[i] = -qm[i][sp]/m*(PoisPot[i+1]-PoisPot[i-1])/(2.0*dx[i]);
  }
  E[nX-1] = -qm[nX-1][sp]/m*(PoisPot[0] - PoisPot[nX-2])/(2.0*dx[nX-1]);  

  for(l=1;l<nX+1;l++)  
    if(E[l] > 0) {
      CFL_NUM = dt*E[l]/h_v;
      //#pragma omp parallel for private(i,j,k,index)
      for(i=0;i<N;i++)
	for(j=0;j<N;j++) 
	  for(k=0;k<N;k++) {	
	    index = k + N*(j + N*i);
	    if(i != 0)
	      f_conv[l][sp][index] = f[l][sp][index] - CFL_NUM*(f[l][sp][index] - f[l][sp][k + N*(j+ N*(i-1))]); 
	    else
	      f_conv[l][sp][index] = f[l][sp][index] - CFL_NUM*(f[l][sp][index]); 
	  }
    }
    else {
      CFL_NUM = dt*E[l]/h_v;
      //#pragma omp parallel for private(i,j,k,index)
      for(i=0;i<N;i++)
	for(j=0;j<N;j++) 
	  for(k=0;k<N;k++) {	
	    index = k + N*(j + N*i);      
	    if(i != (N-1))
	      f_conv[l][sp][index] = f[l][sp][index] - CFL_NUM*( f[l][sp][k + N*(j + N*(i+1))] - f[l][sp][k + N*(j + N*i)]);
	    else
	      f_conv[l][sp][index] = f[l][sp][index] + CFL_NUM*(f[l][sp][index]);
	  }
    }

  free(E);
}


//Computes x portion of second order upwind solution, with minmod
void upwindTwo_x(double ***f, double ***f_conv, double *v, int sp) {
  int i,j,k,l;
  int index;
  double slope[3];
  double CFL_NUM;

  fillGhostCellsPeriodic_secondorder(f,sp);

  //main upwinding
  double f_l, f_r, f_ll, f_rr, x_l, x_r, x_ll, x_rr, dx_l, dx_r;

  for(l=2;l<nX+2;l++) {    
    x_l  = x[l-1];      
    x_ll = x[l-2];
    dx_l = dx[l-1];
    
    x_r  = x[l+1];
    x_rr = x[l+2];
    dx_r = dx[l+1];

    //#pragma omp parallel for private(i,j,k,index,CFL_NUM,f_l,f_ll,f_r,slope)

    for(i=N/2;i<N;i++) {	
      CFL_NUM = dt*v[i]/dx[l];      
      for(j=0;j<N;j++) 
	for(k=0;k<N;k++)  {
	  index = k + N*(j + N*i);
	  //upwind coming from the left
	  //generate the local slopes
	  f_l  = f[l-1][sp][index];
	  f_ll = f[l-2][sp][index];
	  
	  f_r  = f[l+1][sp][index];
	  
	  slope[1] = minmod(2.0*(f[l][sp][index] - f_l            )/(x[l] - x_l ),
			    2.0*(f_r             - f[l][sp][index])/(x_r  - x[l]), 
			        (f_r             - f_l            )/(x_r  - x_l));
	  slope[0] = minmod(2.0*(f_l             - f_ll)/(x_l  - x_ll),
			    2.0*(f[l][sp][index] - f_l )/(x[l] - x_l), 
			        (f[l][sp][index] - f_ll)/(x[l] - x_ll));
	  	  
	  f_conv[l][sp][index] =  - CFL_NUM*(f[l][sp][index] + 0.5*(dx[l] - v[i]*dt)*slope[1] - 
					    (f_l             + 0.5*(dx_l  - v[i]*dt)*slope[0]));
	}
    }

    //#pragma omp parallel for private(i,j,k,index,CFL_NUM,f_l,f_rr,f_r,slope)
    for(i=0;i<N/2;i++) {
      CFL_NUM = dt*v[i]/dx[l];      
      for(j=0;j<N;j++) 
	for(k=0;k<N;k++)  {
	  index = k + N*(j + N*i);
	  //upwind coming from the right
	  //generate the local slopes
	  f_l  = f[l-1][sp][index];

	  f_r  = f[l+1][sp][index];
	  f_rr = f[l+2][sp][index];

	  slope[2] = minmod(2.0*(f_r  - f[l][sp][index])/(x_r - x[l]),
			    2.0*(f_rr - f_r            )/(x_rr - x_r), 
			        (f_rr - f[l][sp][index])/(x_rr - x[l]));

	  slope[1] = minmod(2.0*(f[l][sp][index] - f_l            )/(x[l] - x_l),
			    2.0*(f_r             - f[l][sp][index])/(x_r - x[l]), 
			        (f_r             - f_l            )/(x_r - x_l));
	  
	  //the upwinding
	  f_conv[l][sp][index] = - CFL_NUM*(f_r             - 0.5*(dx_r  + v[i]*dt)*slope[2] - 
					   (f[l][sp][index] - 0.5*(dx[l] + v[i]*dt)*slope[1]));		    	 
	}
    }
  }		
	
}

//Computes seconrd order upwind term for the field acceleration term
void upwindTwo_v(double ***f, double ***f_conv, double *PoisPot, double **qm, double m, double *v, int sp) {
  int i,j,k,l;
  int index;
  double slope[3];
  double CFL_NUM;

  h_v = v[1]-v[0];

  //Set electric field acceleration term
  //This is given by
  // - Z (e_0 E) / m
  //
  // e_0 E = (e phi)_x is given in units of ergs
  // Z is the partial ionization, provided as qm in the input (no unit)
  // m is the mass of the species in grams
  //
  //Thus the acceleration term units are cm / s^2

  double *E = malloc(nX*sizeof(double));
  E[0] = -qm[0][sp]/m*(PoisPot[1]-PoisPot[nX-1])/(2.0*dx[0]);
  for(i=1;i<nX-1;i++) {
    E[i] = -qm[i][sp]/m*(PoisPot[i+1]-PoisPot[i-1])/(2.0*dx[i]);
  }
  E[nX-1] = -qm[nX-1][sp]/m*(PoisPot[0] - PoisPot[nX-2])/(2.0*dx[nX-1]);  

  if (fabs(v[2] - v[1] -  h_v) > 1e-6) {
    printf("Error - you must use a uniform velocity grid for second order advection\n Please set discret=0 in your input file.\n");
    exit(1);
  }

  //main upwinding
  double f_l, f_r, f_ll, f_rr, x_l, x_r, x_ll, x_rr, dx_l, dx_r;
  double up2, up1, neut, down1, down2;

  for(l=0;l<nX;l++) { 
    CFL_NUM = dt*E[l]/h_v;   
    //#pragma omp parallel for private(i,j,k,index,neut,up1,up2,down1,down2,slope)    
    for(i=0;i<N;i++) 
      for(j=0;j<N;j++)
	for(k=0;k<N;k++) {
	  index = k + N*(j + N*i);
	  //Poisson part
	  if(CFL_NUM > 0) {	  
	    neut  = f[l][sp][index];
	    if(i == 0) {
	      up1   = f[l][sp][k + N*(j + N*(i+1))];
	      down1 = 0.0;
	      down2 = 0.0;
	    }
	    else if(i == 1) {
	      up1   = f[l][sp][k + N*(j + N*(i+1))];
	      down1 = f[l][sp][k + N*(j + N*(i-1))];
	      down2 = 0.0;
	    }
	    else if(i == N-1) {
	      up1   = 0.0;
	      down1 = f[l][sp][k + N*(j + N*(i-1))];
	      down2 = f[l][sp][k + N*(j + N*(i-2))];	    
	    }
	    else {
	      up1   = f[l][sp][k + N*(j + N*(i+1))];
	      down1 = f[l][sp][k + N*(j + N*(i-1))];
	      down2 = f[l][sp][k + N*(j + N*(i-2))];	    
	    }
	    
	    slope[1] = minmod(2.0*(up1-neut)    / h_v,
			      2.0*(neut-down1)  / h_v, 
			      (up1-down1)   / (2.0*h_v));
	    slope[0] = minmod(2.0*(neut-down1)  / h_v,
			      2.0*(down1-down2) / h_v,
			      (neut-down2)  / (2.0*h_v));
	    
	    if(i != 0) 
	      f_conv[l][sp][index] = - CFL_NUM*(f[l][sp][index]               + 0.5*(h_v - E[l]*dt)*slope[1] -
					       (f[l][sp][k + N*(j + N*(i-1))] + 0.5*(h_v - E[l]*dt)*slope[0]));
	    else
	      f_conv[l][sp][index] = - CFL_NUM*(f[l][sp][index] + 0.5*(h_v - E[l]*dt)*slope[1]);	    
	  }
	  else {
	    neut  = f[l][sp][index];
	    if(i == 0) {
	      up2   = f[l][sp][k + N*(j + N*(i+2))];	    
	      up1   = f[l][sp][k + N*(j + N*(i+1))];
	      down1 = 0.0;
	    }
	    else if(i == N-1) {
	      up2   = 0.0;	    
	      up1   = 0.0;
	      down1 = f[l][sp][k + N*(j + N*(i-1))];
	    }
	    else if(i == N-2) {
	      up2   = 0.0;	    
	      up1   = f[l][sp][k + N*(j + N*(i+1))];
	      down1 = f[l][sp][k + N*(j + N*(i-1))];
	    }
	    else {
	      up2   = f[l][sp][k + N*(j + N*(i+2))];
	      up1   = f[l][sp][k + N*(j + N*(i+1))];
	      down1 = f[l][sp][k + N*(j + N*(i-1))];
	    }
	    
	    
	    slope[2] = minmod(2.0*(up2-up1)    / h_v,
			      2.0*(up1-neut)   / h_v, 
			      (up2-neut)   / (2.0*h_v));
	    slope[1] = minmod(2.0*(up1-neut)   / h_v,
			      2.0*(neut-down1) / h_v,
			      (up1-down1)  / (2.0*h_v));
	    
	    if(i != (N-1))
	      f_conv[l][sp][index] = - CFL_NUM*(f[l][sp][k + N*(j + N*(i+1))] - 0.5*(h_v + E[l]*dt)*slope[2] -
					       (f[l][sp][index]               - 0.5*(h_v + E[l]*dt)*slope[1]));	  
	    else
	      f_conv[l][sp][index] =   CFL_NUM*(f[l][sp][index] - 0.5*(h_v + E[l]*dt)*slope[1]);	  
	  }
	  
	}
  }     

  free(E);  
}  



//does a first order splitting, solving v first and returns the updated f
void advectOne(double ***f, double *PoisPot, double **qm, double m, int sp) {
  upwindOne_v(f, f_star, PoisPot, qm, m, c[sp],sp);
  upwindOne_x(f_star, f, c[sp],sp);
}

//Defunct, do not use. This is a placeholder for a first order time splitting, but with high-resolution methods for each piece
void advectTwo(double ***f, double *PoisPot, double **qm, double m, int sp) {
  upwindTwo_v(f, f_star, PoisPot, qm, m, c[sp], sp);   
  upwindTwo_x(f_star, f, c[sp], sp);   
}

//This returns the *update* to f as f_conv for the x transport piece
//used in second order time stepping method
void advectTwo_x(double ***f, double ***f_conv, int sp) {
  upwindTwo_x(f, f_conv, c[sp], sp);   
}

//This returns the *update* to f as f_conv for the v transport piece
//used in second order time stepping method
void advectTwo_v(double ***f, double ***f_conv, double *PoisPot, double **qm, double m, int sp) {
  upwindTwo_v(f, f_conv, PoisPot, qm, m, c[sp], sp);   
}


void dealloc_trans() {
  int i,l;
  for(l=0;l<nX;l++) {
    for(i=0;i<ns;i++) 
      free(f_star[l][i]);
    free(f_star[l]);
  }
  free(f_star);
}
