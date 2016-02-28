#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "riemann.h"

/* Punto 1:
Escriba un código en C que resuelva el problema del shock-tube con
un esquema de Godunov de primer orden. Utilice las mismas variables que
se inicializaron en el [ejemplo con el esquema Lax-Friedrichs]
(https://github.com/ComputoCienciasUniandes/MetodosComputacionalesAvanzados/blob/master/weeks/03/code/shocktube.c).

C adaptation from C++ code written by Richard J. Gonsalves.
Based on Godunov's explanation here: http://www.physics.buffalo.edu/phy411-506/topic7/topic7-lec2.pdf
(Specially pages 1-9)
*/

double L = 4.0;                   // length of shock tube
double gama = 1.4;                // ratio of specific heats
int N = 5000;                   // number of grid points

double CFL = 0.4;                 // Courant-Friedrichs-Lewy number

double **U = NULL;                // solution with 3 components
double **F = NULL;                // flux with 3 components

double h;                         // lattice spacing
double tau;                       // time step
double c;                         // speed

int step;

void allocate();
double cMax();
void initialize();
void boundaryConditions(double **U);
void upwindGodunovStep();
void solve(double tMax);

int main()
{
  solve(1.0);

  return 0;
}

void allocate() {

  int j;

  U = malloc(N * sizeof(double *));
  F = malloc(N * sizeof(double *));

  for (j = 0; j < N; j++) {
    U[j] = malloc(3 * sizeof(double));
    F[j] = malloc(3 * sizeof(double));
  }
}

double cMax() {

  double uMax = 0;
  double rho, u, p, c;
  int i;

  for (i = 0; i < N; i++) {
    if (U[i][0] == 0)
	   continue;

    rho = U[i][0];
    u = U[i][1]/rho;
    p = (U[i][2]-rho*u*u/2)*(gama-1);

    c = sqrt(gama*fabs(p)/rho);

    if (uMax < (c + fabs(u)))
	   uMax = c + fabs(u);
  }

  return uMax;
}

void initialize() {

  int j;
  double rho,p,u,e;
  allocate();

  h = 1.0 * L / (N - 1);

  for (j = 0; j < N; j++) {

    rho = 1;
    p = 1;
    u = 0;

    if (j > N / 2){
      rho = 0.125;
      p = 0.1;
    }

    e = p/(gama-1) + rho*u*u/2.0;

    U[j][0] = rho;
    U[j][1] = rho*u;
    U[j][2] = e;

    F[j][0] = rho*u;
    F[j][1] = rho*u*u+p;
    F[j][2] = u*(e+p);
  }

  tau = CFL*h/cMax();
  step = 0;
}

void boundaryConditions(double **U) {
  // reflection boundary conditions at the tube ends
  U[0][0] = U[1][0];
  U[0][1] = -U[1][1];
  U[0][2] = U[1][2];

  U[N-1][0] = U[N-2][0];
  U[N-1][1] = -U[N-2][1];
  U[N-1][2] = U[N-2][2];
}

void upwindGodunovStep() {

  int i, j;

  // find fluxes using Riemann solver
  for (j = 0; j < N - 1; j++){
    Riemann(U[j], U[j + 1], F[j]);
  }

  // update U
  for (j = 1; j < N - 1; j++){
    for (i = 0; i < 3; i++){
      U[j][i] -= tau/h*(F[j][i] - F[j-1][i]);
    }
  }
}

void solve(double tMax)
{
  initialize();

  double t = 0.0;

  int j;

  double rho, u, e, P;

  tau = CFL*h/cMax();

  while (t < tMax) {

  	boundaryConditions(U);

  	tau = CFL*h/cMax();

    upwindGodunovStep();

  	t += tau;

    // rho = U[98][0];
    // u = U[98][1]/rho;
    // e = U[98][2];
    // P = (gama-1.0)*(e-rho*u*u/2.0);
    // printf("%f\t%.20f\t%.20f\t%.20f\t%.20f\n", t, rho, u, e, P);
    // rho = U[102][0];
    // u = U[102][1]/rho;
    // e = U[102][2];
    // P = (gama-1.0)*(e-rho*u*u/2.0);
    // printf("%f\t%.20f\t%.20f\t%.20f\t%.20f\n", t, rho, u, e, P);
  }

  //File to plot last step
  FILE *final_step_file;
  final_step_file = fopen("final_step.dat", "w");
  double current_x;

  for (j = 0; j < N; j++) {

    rho = U[j][0];
    u = U[j][1]/rho;
    e = U[j][2];
    P = (gama-1.0)*(e-rho*u*u/2.0);

    current_x = j*h;

    fprintf(final_step_file, "%f\t%.20f\t%.20f\t%.20f\t%.20f\n", current_x, rho, u, e, P);
  }
}
