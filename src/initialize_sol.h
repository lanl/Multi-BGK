
void initialize_sol_inhom(double ***f, int nint, double *int_loc, double *rho_in, double *v_in, double *T_in, int Nx, double *x, int nspec, int nV, int order, double **c, double *m, double **n_oned, double ***v_oned, double **T_oned);

void initialize_sol_load_inhom_file(int Nx, int nspec, double **n_oned, double ***v_oned, double **T_oned, char *input_file_data_filename);

void initialize_sol_inhom_file(double ***f, int Nx, int nspec, int nV, int order, double **c, double *m, double **n_oned, double ***v_oned, double **T_oned);

