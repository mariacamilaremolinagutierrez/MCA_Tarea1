/*
  Translated from a C++ Riemann solver written by Richard
  J. Gonsalves,  which in turn is rewrite of a Riemann solver
  found in Laney's textbook Computational Gasdynamics (Cambridge
  University Press.)
 */

#ifndef RIEMANN_H_INCLUDED
#define RIEMANN_H_INCLUDED

void Riemann(double *U4, double *U1, double *F) {

  const double gamma = 1.4;
  const double g1 = (gamma-1) / (2*gamma);
  const double g2 = (gamma+1) / (2*gamma);
  const double g3 = (gamma+1) / (gamma-1);
  const double tol = 1e-10;

  int i;
  double x, y, z, fz;
  double a1, a2, a3, a4;
  double p1, p2, p3, p4;
  double u1, u2, u3, u4;
  double s1, s2, s3, s4;
  double rho1, rho2, rho3, rho4;

  // compute primitive variables
  rho1 = U1[0];
  u1 = U1[1]/rho1;
  p1 = (U1[2]-rho1*u1*u1/2)*(gamma-1);
  a1 = sqrt(gamma*p1/rho1);

  rho4 = U4[0];
  u4 = U4[1]/rho4;
  p4 = (U4[2]-rho4*u4*u4/2)*(gamma-1);
  a4 = sqrt(gamma*p4/rho4);

  // apply the bisection method to find p2
  x = 0.001; //x=p2/p1
  y = 1.2*p4/p1; //y=p2/p1

  while(y-x > tol){

    z = x+(y-x)/2.0;

    fz = pow(p4/p1,1.0/g1) - pow(z, 1.0/g1)*(1 + ((gamma-1)/(2*a4))*(u4-u1-(a1/gamma)*(z-1)/sqrt(g2*(z-1)+1)));

    if (fz >0) {
      x=z;
    }
    else{
      y=z;
    }
  }

  // Compute shock
  p2 = p1 * x;
  u2 = u1 + (a1/gamma)*((x-1)/sqrt(g2*(x-1)+1));
  a2 = a1 * sqrt(x*(g3+x)/(1+g3*x));
  rho2 = gamma*p2/(a2*a2);

  //s = u1 + a1*sqrt(g2*(x-1)+1);
  s1 = u1 + a1*sqrt(g2*(x-1)+1);

  // Compute contact
  p3 = p2;
  u3 = u2;
  a3 = (u4+2.0*a4/(gamma-1)-u3)*(gamma-1)/2.0;
  rho3 = gamma*p3/(a3*a3);

  s2 = u2;

  // Compute expansion
  s3 = u3 - a3;
  s4 = u4 - a4;

  // Compute fluxes
  double f1, f2, f3;
  double a, u, p, rho;

  if(s4 > 0) {
    f1 = rho4*u4;
    f2 = rho4*u4*u4 + p4;
    f3 = 0.5*rho4*u4*u4*u4 + rho4*a4*a4*u4/(gamma-1.);
  }

  else if (s3 > 0) {
    u = ((gamma-1.0)*u4+2.0*a4)/(gamma+1.);
    a = u;
    p = p4*pow(a/a4, 2.0*gamma/(gamma-1.));

    if (a < 0 || p < 0) {
      printf("Negative a or p in Riemann");
    }

    rho = gamma*p/(a*a);
    f1 = rho*u;
    f2 = rho*u*u + p;
    f3 = 0.5*rho*u*u*u + rho*a*a*u/(gamma-1.0);
  }

  else if (s2 > 0) {
    f1 = rho3*u3;
    f2 = rho3*u3*u3 + p3;
    f3 =  0.5*rho3*u3*u3*u3 + rho3*a3*a3*u3/(gamma-1.0);
  }

  else if (s1 > 0) {
    f1 = rho2*u2;
    f2 = rho2*u2*u2 + p2;
    f3 = 0.5*rho2*u2*u2*u2 + rho2*a2*a2*u2/(gamma-1.0);
  }

  else {
    f1 = rho1*u1;
    f2 = rho1*u1*u1 + p1;
    f3 = 0.5*rho1*u1*u1*u1 + rho1*a1*a1*u1/(gamma-1.0);
  }

  F[0] = f1;
  F[1] = f2;
  F[2] = f3;
}

#endif /* RIEMANN_H_INCLUDED */
