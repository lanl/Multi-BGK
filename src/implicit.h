void implicitGetVelocitiesTemperaturesLinear(double *n, double **v, double *T,
                                             double **nu, double *m, double dt,
                                             int nspec, int conserveFlag,
                                             double **vnew, double ***vmix_new,
                                             double *Tnew, double **Tmix_new);

void implicitGetVelocitiesTemperaturesNonlinear(
    double *n, double **v, double *T, double *m, double dt, int nspec,
    int conserveFlag, double **vnew, double ***vmix_new, double *Tnew,
    double **Tmix_new);

void implicitUpdateDistributionsLinear(double **f, double *n, double **v,
                                       double *T, double **nu, double *m,
                                       double dt, int nspec, int conserveFlag,
                                       int Nv, double **fnew);

void implicitUpdateDistributionsNonlinear(double **f, double *n, double **v,
                                          double *T, double **nu, double *m,
                                          double dt, int nspec,
                                          int conserveFlag, int Nv,
                                          double **fnew);
