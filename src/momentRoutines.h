void initialize_moments(int nodes, int nspec, double **vel, double **wts);

double getDensity(double *in, int sp);

void getBulkVel(double *in, double *out, double n, int sp);

double getTemp(double m, double n, double *u, double *in, int sp);

double getDirTemp(double m, double n, double *u, double *in, int sp, int dir);

double getH(double n, double *in, int sp);

void getMarginal(double *in, double *marginal, int sp);

//computes n Theta, the pressure/stress tensor
void getStressTensor(double m, double n, double *u, double *in, double *out, int sp);
