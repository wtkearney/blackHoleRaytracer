
#ifndef ____rk4__
#define ____rk4__

double rk4 ( double t0, double u0, double dt, double f ( double t, double u ) );
double *rk4vec ( double t0, int n, double u0[], double dt, 
  double *f ( double t, int n, double u[] ) );
void timestamp ( void );

#endif /* defined(____rk4__) */