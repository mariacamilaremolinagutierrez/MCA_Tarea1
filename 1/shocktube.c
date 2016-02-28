#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
int N = 5000;                     // number of grid points

double CFL = 0.4;                 // Courant-Friedrichs-Lewy number
double nu = 0.0;                  // artificial viscosity coefficient

double **U = NULL;                // solution with 3 components

double h;                         // lattice spacing
double tau;                       // time step
double c;                         // speed

int step;

void allocate();
double cMax();
void initialize();
void boundaryConditions(double **U);
void upwindGodunovStep();
void solve(double tMax, char *filename, int plots);

int main()
{
  solve(1.0, "UpwindGodunov", 5);

  return 0;
}

void allocate() {

  int j;

  U = malloc(N * sizeof(double *));

  for (j = 0; j < N; j++) {
    U[j] = malloc(3 * sizeof(double));
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
    U[j][1] = rho * u;
    U[j][2] = e;
  }

  c = cMax();

  tau = CFL*h/c;
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

  for (j = 1; j < N - 1; j++){
    for (i = 0; i < 3; i++){
      if(c>0){
        U[j][i] -= CFL*tau/h*(U[j][i] - U[j-1][i]);
      }
      else{
        U[j][i] -= CFL*tau/h*(U[j+1][i] - U[j][i]);
      }
    }
  }
}

void solve(double tMax, char *filename, int plots)
{
  initialize();

  double t = 0.0;
  int step = 0;
  int plot = 0;

  FILE *out;
  char filename_tmp[1024];

  int j;

  double rho, u, e, P;

  tau = CFL*h/c;

  while(plot<=plots) {

    sprintf(filename_tmp, "%s/%s_step_%d.dat", filename, filename, plot);

    if(!(out = fopen(filename_tmp, "w"))){
    	fprintf(stderr, "problem opening file %s\n", filename);
    	exit(1);
    }

    // write solution in plot files and print
    double rho_avg = 0.0, u_avg = 0.0, e_avg = 0.0, P_avg = 0.0;

    for (j = 0; j < N; j++) {
    	rho = U[j][0];
    	u = U[j][1] / U[j][0];
    	e = U[j][2];
    	P = (U[j][2] - U[j][1] * U[j][1] / U[j][0] / 2) * (gama - 1.0);

    	rho_avg += rho;
    	u_avg += u;
    	e_avg += e;
    	P_avg += P;

    	fprintf(out, "%d\t%.20f\t%.20f\t%.20f\t%.20f\n", j, rho, u, e, P);
    }

    fclose(out);

    if (rho_avg != 0.0) rho_avg /= N;
    if (u_avg != 0.0)   u_avg /= N;
    if (e_avg != 0.0)   e_avg /= N;
    if (P_avg != 0.0)   P_avg /= N;

    fprintf(stdout,"Step %d Time %f\tRho_avg %f\t u_avg %f\t e_avg %f\t P_avg %f\n", step, t,rho_avg,u_avg,e_avg,P_avg);

    plot++;

    while (t < tMax*plot / (double)(plots)) {
    	boundaryConditions(U);
    	tau = CFL*h/c;
      upwindGodunovStep();
    	t += tau;
    	step++;
    }
  }

  //File to plot las step
  FILE *final_step_file;
  final_step_file = fopen("final_step.dat", "w");
  double current_x;

  for (j = 0; j < N; j++) {
    rho = U[j][0];
    u = U[j][1] / U[j][0];
    e = U[j][2];
    P = (U[j][2] - U[j][1] * U[j][1] / U[j][0] / 2) * (gama - 1.0);

    current_x = j*h;

    fprintf(final_step_file, "%f\t%.20f\t%.20f\t%.20f\t%.20f\n", current_x, rho, u, e, P);
  }
}
