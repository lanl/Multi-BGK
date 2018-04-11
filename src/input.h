#ifndef _INPUT_H
#define _INPUT_H

void read_input(int *nspec, 
		int *dims, 
		int *Nx, 		
		double *Lx, 
		int *Nv, 
		double *v_sigma,
		int *discret, 
		int *poissFlavor, 
		double **m, 
		double **Z, 
		int *order, 
		int *im_ex, 
		double *dt, 
		double *tfinal, 
		int *numint,
		double **intervalLimits,
		double **ndens_in,
		double **velo_in,
		double **T_in,
		int *ecouple, 
		double *Te_start,
		double *Te_end,
		int *CL_type, 
		int *ion_type, 
		int *MT_or_TR,  
		double **n, 
		double **u, 
		double **T, 
		int *dataFreq, 
		int *outputDist,
		double *RHS_tol,
		int *BGK_type, 
		double *beta, 
		int *hydro_flag,
		int *TNBFlag,
		char *inputFilename);

void set_default_values(int *Nx, 
			double *Lx, 
			int *Nv,
			double *v_sigma,
			int *order, 
			int *discret, 
			int *im_ex, 
			int *poissFlavor, 
			int *ecouple, 
			double *Te_start,
			double *Te_end,
			int *CL_type, 
			int *ion_type, 
			int *MT_or_TR, 
			double *dt, 
			double *tfinal, 
			int *BGK_type, 
			double *beta, 
			int *hydro_flag,
			int *TNBFlag,
			int *dataFreq,
			int *dumpDist,
			double *RHS_tol);


void check_input(const int *flag);

void read_line(FILE *file, char line[80]);
 
size_t read_int(FILE *file);

double read_double(FILE *file);
 
void read_line_no_adv(FILE *file, char line[80]);

size_t read_int_no_adv(FILE *file);

double read_double_no_adv(FILE *file);

#endif
