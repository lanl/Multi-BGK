void initialize_moments(int nodes, int nspec, double **vel, double **wts);

double getDensity(double *in, int sp);

void getBulkVel(double *in, double *out, double n, int sp);

double getTemp(double m, double n, double *u, double *in, int sp);

double getH(double n, double *in, int sp);

void getSecondMoment_Tensor(double m, double *in, double **out, int sp);

void getThirdMoment_Vector(double m, double *in, double *out, int sp);
