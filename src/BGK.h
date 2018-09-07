void GetMaxwell(double m, double n, double *v, double T, double *M, int sp);

void initialize_BGK(double ns, int numV, double *mass, double **vels, int order,
                    int ec, int CL, int itype, int MorT, int tauFlag,
                    char *filename);

void BGK_ex(double **f, double **f_out, double *Z, double dt, double Te);

void BGK_im(double **f, double *Z, double dt, double Te);

void BGK_Greene(double **f, double **f_out, double *Z, double dt, double beta,
                double Te);

void BGK_NRL(double **f, double **f_out, double *Z, double dt, double Te);

void BGK_norm(double **f, double **f_err, double *Z, double dt, double Te);
