void initialize_fluid_electron(int Nx_rank_in, int nspec_in, double dt_in);

void fluid_electron_coll_only(double **n, double *Te, double **Ti);

void fluid_electron_get_eE(int Nx, double *n_e, double dx, double *eE, double *Te);

void fluid_electron_ohmic_term(double **n, double ***v_e, double *T_e, double *eE, double **Z);
