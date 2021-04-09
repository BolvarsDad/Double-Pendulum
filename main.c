// C program to numerically solve ODEs of a double pendulum integrated with 4th order Runge Kutta.
// README for instructions.
// Compilation example: ./a.exe 0.0 10.0 90.0 0.00 -10.0 0.0 1000 > outfile.txt.
// Where given arguments represent initial values of the pendulum.
// TMIN TMAX TH10 W10 TH20 W20 NSTEP.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Hardwired parameters.

#define N 4	        // Numer of equations to be solved
#define G 9.82	    // Gravity 	             	 (m/s^2)
#define L1 1.0	    // Length of first rod 	 (m)
#define L2 1.0	    // Length of second rod  	 (m)
#define M1 1.0	    // Mass of first bob 	 (kg)
#define M2 1.0	    // Mass of second bob 	 (kg)

void 
rungeKutta(float xin, float yin[], float yout[], float h);

void 
derivs(float xin, float yin[], float dydx[]);

// Main function

int 
main(int argc, char **argv)
{
  int i = 0, NSTEP;
  float h, TMIN, TMAX, TH10, W10, TH20, W20;
  float yin[N], yout[N];
  float *t, *th1, *th2, *w1, *w2;

  // obtain command line values
  // atof takes a string argument and converts it to a floating point number (type double)
  // atoi takes a string argument and converts it to an integer (type int)

  TMIN  = atof(argv[1]);	    // Starting time 	             (s).
  TMAX  = atof(argv[2]);	    // Ending time   	             (s).
  TH10  = atof(argv[3]);	    // Initial angle 1 	             (d).
  W10   = atof(argv[4]);	    // Initial angular velocity 1    (d/s).
  TH20  = atof(argv[5]);	    // Initial angle 2	             (d).
  W20   = atof(argv[6]);	    // Initial angular velocity 2    (d/s).
  NSTEP = atoi(argv[7]);	    // Number of integrations per step.

  // allocate memory for arrays of values of time, angles 1 and 2,
  // and angular velocities 1 and 2 respectively 

  t   = (float *) malloc(NSTEP * sizeof(float)); 
  th1 = (float *) malloc(NSTEP * sizeof(float)); 
  w1  = (float *) malloc(NSTEP * sizeof(float));
  th2 = (float *) malloc(NSTEP * sizeof(float));
  w2  = (float *) malloc(NSTEP * sizeof(float));

  // stepsize for integration

  h = (TMAX - TMIN) / (NSTEP - 1.0);

  // Define array of t values

  for (i = 0; i < NSTEP; i++)
  {
	  t[i] = TMIN + h * i;
  }

  // initial values - convert all angles to radians
  // (r) = radians
  
  th1[0]  = TH10 * M_PI / 180.0;	      // Theta 1 	        (r)
   w1[0]  =  W10 * M_PI / 180.0;	      // Theta velocity 1 	(r/s)
  th2[0]  = TH20 * M_PI / 180.0;	      // Theta 2 	        (r)
   w2[0]  =  W20 * M_PI / 180.0; 	      // Theta velocity 2 	(r/s)

  // perform the integration

  printf("%f %f %f %f %f\n", t[0], th1[0], w1[0], th2[0], w2[0]);

  // for loop to print integrated results

  for (i = 0; i < NSTEP - 1; i++)
  { 
    yin[0] = th1[i];
    yin[1] =  w1[i];
    yin[2] = th2[i];
    yin[3] =  w2[i];

    rungeKutta(t[i], yin, yout, h);

    th1[i+1]  = yout[0];
     w1[i+1]  = yout[1];
    th2[i+1]  = yout[2];
     w2[i+1]  = yout[3];

    printf("%f %f %f %f %f\n", t[i+1], th1[i+1], w1[i+1], th2[i+1], w2[i+1]);
  }

  return 0;
}

void
derivs(float xin, float yin[], float dydx[])
{

  // function to fill array of derivatives dydx at xin

  float den1, den2, del;

  dydx[0] = yin[1]; 
  
  del     = yin[2]-yin[0];
  den1    = (M1 + M2) * L1 - M2 * L1 * cos(del) * cos(del);
  dydx[1] = (M2 * L1 * yin[1] * yin[1] * sin(del) * cos(del)
    	    + M2 * G * sin(yin[2]) * cos(del) + M2 * L2 * yin[3] * yin[3] * sin(del)
    	    - (M1+M2) * G * sin(yin[0]))/den1;

  dydx[2] = yin[3];

  den2    = (L2 / L1) * den1;
  dydx[3] = (-M2 * L2 * yin[3] * yin[3] * sin(del) * cos(del)
    	    + (M1 + M2) * G * sin(yin[0]) * cos(del) 
    	    - (M1 + M2) * L1 * yin[1] * yin[1] * sin(del)
    	    - (M1 + M2) * G * sin(yin[2])) / den2;

  return;
}

void
rungeKutta(float xin, float yin[], float yout[], float h)
{

  // fourth order Runge-Kutta
 
  int i;
  float hh, xh, dydx[N], dydxt[N], yt[N], k1[N], k2[N], k3[N], k4[N];
  
  hh = 0.5 * h;
  xh = xin + hh;
  
  // First step

  derivs(xin, yin, dydx);
  for (i = 0; i < N; i++) 
  {
    k1[i] = h * dydx[i];
    yt[i] = yin[i] + 0.5 * k1[i];
  }

  // Second Step

  derivs(xh, yt, dydxt);
  for (i = 0; i < N; i++)
  {
    k2[i] = h * dydxt[i];
    yt[i] = yin[i] + 0.5 * k2[i];
  }   

  // Third step

  derivs(xh, yt, dydxt);
  for (i = 0; i < N; i++)
  {
    k3[i] = h * dydxt[i];
    yt[i] = yin[i] + k3[i];
  }

  // Fourth step

  derivs(xin + h, yt, dydxt);
  for (i = 0; i < N; i++)
  {
    k4[i]   = h * dydxt[i];
    yout[i] = yin[i] + k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0;
  }

  return;
}
