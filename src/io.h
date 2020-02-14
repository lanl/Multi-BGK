
void io_init_homog(int numV, int numS, double **velos);

void io_init_inhomog(int numX, int numV, int numS, double **velos);

void store_distributions_homog(double **f, double t, int nT, char *fileName);

void load_distributions_homog(double **f, char *fileName);

void store_distributions_inhomog(double ***f, char *fileName, int t);

void store_grid(char *fileName);

void load_grid_restart(double *Lv, double *t, int *nT, char *fileName);

void load_taus_homog(double **nu, char *fileName);

void load_diffusion_homog(double **nu, char *fileName);

#ifdef ALDR_ON

void request_aldr_single(double *n, double *T, double *Z, char *tag, char *dbfile, double **D_ij);

void request_aldr_batch(double **n, double **T, double **Z, char *tag, char *dbfile, double ***D_ij, int *provenance);

void test_aldr();

void io_init_db(char *filename);

void io_close_db();

#endif
