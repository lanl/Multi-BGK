
void GetMaxwell(double m, double n, double *v, double T, double *M, int sp);

void initialize_BGK(double ns, int numV, double *mass, double **vels, int order, int ec, int CL, int itype, int MorT, int tauFlag, char *filename);

void BGK_ex(double **f, double **f_out, double *Z, double dt, double Te);

void BGK_im_linear(double **f, double **f_out, double *Z, double dt, double Te);

void BGK_im_nonlinear(double **f, double **f_out, double *Z, double dt, double Te);

void BGK_Greene(double **f, double **f_out, double *Z, double dt, double beta, double Te);

void BGK_NRL(double **f, double **f_out, double *Z, double dt, double Te);

void BGK_norm(double **f, double **f_err, double *Z, double dt, double Te);

#ifdef ALDR_ON
void get_diffusion_from_MD_0d(double *n, double *T, double *Z, char *tag, char *dbname);

void set_diffusion_from_MD_1d(double **Dij_in);

#endif
