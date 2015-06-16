# include <stdlib.h>
# include <stdio.h>
# include <time.h>

#include "graphics.h"

/******************************************************************************/

double rk4 ( double t0, double u0, double dt, double f ( double t, double u ) )

/******************************************************************************/
/*
  Purpose:
 
    RK4 takes one Runge-Kutta step for a scalar ODE.

  Discussion:

    It is assumed that an initial value problem, of the form

      du/dt = f ( t, u )
      u(t0) = u0

    is being solved.

    If the user can supply current values of t, u, a stepsize dt, and a
    function to evaluate the derivative, this function can compute the
    fourth-order Runge Kutta estimate to the solution at time t+dt.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, double T0, the current time.

    Input, double U0, the solution estimate at the current time.

    Input, double DT, the time step.

    Input, double F ( double T, double U ), a function which evaluates
    the derivative, or right hand side of the problem.

    Output, double RK4, the fourth-order Runge-Kutta solution estimate
    at time T0+DT.
*/
{
  double f0;
  double f1;
  double f2;
  double f3;
  double t1;
  double t2;
  double t3;
  double u;
  double u1;
  double u2;
  double u3;
/*
  Get four sample values of the derivative.
*/
  f0 = f ( t0, u0 );

  t1 = t0 + dt / 2.0;
  u1 = u0 + dt * f0 / 2.0;
  f1 = f ( t1, u1 );

  t2 = t0 + dt / 2.0;
  u2 = u0 + dt * f1 / 2.0;
  f2 = f ( t2, u2 );

  t3 = t0 + dt;
  u3 = u0 + dt * f2;
  f3 = f ( t3, u3 );
/*
  Combine them to estimate the solution.
*/
  u = u0 + dt * ( f0 + 2.0 * f1 + 2.0 * f2 + f3 ) / 6.0;

  return u;
}
/******************************************************************************/

double *rk4vec ( double t0, int m, double u0[], double dt, 
  double *f ( double t, int m, double u[] ) )

/******************************************************************************/
/*
  Purpose:
 
    RK4VEC takes one Runge-Kutta step for a vector ODE.

  Discussion:

    It is assumed that an initial value problem, of the form

      du/dt = f ( t, u )
      u(t0) = u0

    is being solved.

    If the user can supply current values of t, u, a stepsize dt, and a
    function to evaluate the derivative, this function can compute the
    fourth-order Runge Kutta estimate to the solution at time t+dt.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, double T0, the current time.

    Input, int M, the spatial dimension.

    Input, double U0[M], the solution estimate at the current time.

    Input, double DT, the time step.

    Input, double *F ( double T, int M, double U[] ), a function which evaluates
    the derivative, or right hand side of the problem.

    Output, double RK4VEC[M], the fourth-order Runge-Kutta solution estimate
    at time T0+DT.
*/
{
  double *f0;
  double *f1;
  double *f2;
  double *f3;
  int i;
  double t1;
  double t2;
  double t3;
  double *u;
  double *u1;
  double *u2;
  double *u3;
/*
  Get four sample values of the derivative.
*/
  f0 = f ( t0, m, u0 );

  t1 = t0 + dt / 2.0;
  u1 = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    u1[i] = u0[i] + dt * f0[i] / 2.0;
  }
  f1 = f ( t1, m, u1 );

  t2 = t0 + dt / 2.0;
  u2 = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    u2[i] = u0[i] + dt * f1[i] / 2.0;
  }
  f2 = f ( t2, m, u2 );

  t3 = t0 + dt;
  u3 = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
     u3[i] = u0[i] + dt * f2[i];
  }
  f3 = f ( t3, m, u3 );
/*
  Combine them to estimate the solution.
*/
  u = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
     u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;
  }
/*
  Free memory.
*/
  free ( f0 );
  free ( f1 );
  free ( f2 );
  free ( f3 );
  free ( u1 );
  free ( u2 );
  free ( u3 );

  return u;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}


