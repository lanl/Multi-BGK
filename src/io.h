
void io_init_homog(int numV, int numS, double **velos);

void io_init_inhomog(int numX, int numV, int numS, double **velos);

void store_distributions_homog(double **f, double t, int nT, char *fileName);

void load_distributions_homog(double **f, char *fileName);

void store_distributions_inhomog(double ***f, char *fileName, int t);

void store_grid(char *fileName);

void load_grid_restart(double *Lv, double *t, char *fileName);

void load_taus_homog(double **nu, char *fileName);
