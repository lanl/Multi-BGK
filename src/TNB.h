void initializeTNB(int Nv_in, double **c_in, double **wts_in);

double GetTNB_dt(double mu, double *in, double *c1, int sp, int sp2);
double GetReactivity_dt(double mu, double *in,double *in2, int sp, int sp2);

double GetTNB_dd_He(double mu, double *in, double *c1, int sp, int sp2);
double GetReactivity_dd_He(double mu, double *in,double *in2, int sp, int sp2);

double GetTNB_dd_T(double mu, double *in, double *c1, int sp, int sp2);
double GetReactivity_dd_T(double mu, double *in,double *in2, int sp, int sp2);

double GetTNB_tt(double mu, double *in, double *c1, int sp, int sp2);
double GetReactivity_tt(double mu, double *in,double *in2, int sp, int sp2);

