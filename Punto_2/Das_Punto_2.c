#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//Definicion de constantes
//Opciones para ejecutar el codigo
#define OPT_SIMPLETIC 1
#define OPT_RK4 2

//Metodos para el manejo de arreglos y matrices
double* crearArregloCero(int n_points);
double* crearArregloEquiEspaciado(double x_ini,double x_fin,int n_points);
double* copiarArreglo(double* original, int n_points);
double* copiarMatriz(double* original, int n, int m);
double* crear_matriz(int n, int m);

//Metodos para ambos metodos
double hamiltoniano_1(double q_1, double p_1, double epsilon);
double hamiltoniano_3(double q_3, double p_3, double epsilon);
double hamiltoniano_t(double q_1, double p_1,double q_3, double p_3, double epsilon);
double p_punto_3(double q_1,double q_3, double epsilon);
double p_punto_1(double q_1,  double q_3, double epsilon);
double square(double x);
double q_punto(double p);

//Metodos para RK4
void RK4(double big_T, double delta_t, double epsilon, double a, double q_3_0, double p_3_0);
double* RK4_step_1(double delta_t, double q_1, double p_1, double q_3, double p_3, double epsilon);
double* RK4_step_3(double delta_t, double q_1, double p_1, double q_3, double p_3, double epsilon);
void escribirArreglos_RK4(double*t, double*q_1, double*q_3, double*p_1, double*p_3, double* energia, int n_step);

//Metodos para la integracion simpletica
double* gamma_A_3 (double alpha, double q_1, double q_3, double p, double epsilon);
double* gamma_A_1 (double alpha, double q_1, double q_3, double p, double epsilon);
double* gamma_U_3 (double alpha, double q_1, double q_3, double p, double epsilon);
double* gamma_U_1 (double alpha, double q_1, double q_3, double p, double epsilon);
double* labmda_mayuscula_1(double alpha,double q_1, double q_3, double p_1, double epsilon);
double* labmda_mayuscula_3(double alpha,double q_1, double q_3, double p_3, double epsilon);
double* Simpletic_Integration_step(double alpha_1, double alpha_0, double q_1, double p_1, double q_3, double p_3, double epsilon);
void escribirArreglos_simpletico(double*t, double*q_1, double*q_3, double*p_1, double*p_3, double* energia, int n_step);
void Simpletic_Integration(double big_T, double delta_t, double epsilon, double a, double q_3_0, double p_3_0);

/*
 * Metodo que se ejecuta con el programa
 * Se espera que se ejecute con ./Das_Executor.x OPTION_EXECUTOR(1 o 2) T_Final(double) delta_t(double a(double) IC_masa_3(archivo_cond_inic)
 */
int main(int argc, char **argv){
  int opcion;
  opcion = atoi(argv[1]);
  double big_T, delta_t, epsilon, a;
  big_T=atof(argv[2]);
  delta_t=atof(argv[3]);
  
  epsilon = 1.0;
  a = atof(argv[4]);

  //Archivo Condiciones
  int n_condiciones;
  FILE* IC_masa_3;
  if(!(IC_masa_3=fopen(argv[5], "r"))){
    printf("Problem opening file %s\n", argv[5]);
    exit(1);
  }
  //La primera linea del archivo indica el numero de condiciones que se van a utilizar.
  fscanf(IC_masa_3, "%d \n", &n_condiciones);
  int i;
  for(i=0;i<n_condiciones;i++){
    double q_3_0, p_3_0;
    fscanf(IC_masa_3, "%lf %lf\n", &q_3_0, &p_3_0);
    //Mira si se va a ejecutar el Simpletico o el RK4
    if(opcion == OPT_SIMPLETIC){
      Simpletic_Integration(big_T, delta_t, epsilon, a, q_3_0, p_3_0);
    }
    else if (opcion == OPT_RK4){
      RK4(big_T, delta_t, epsilon, a, q_3_0, p_3_0);
    }
  }
  
  
  
  
  return 0;
}

/*
 *---------------------------------------------------------------------------------------------
 * RK4
 */
void RK4(double big_T, double delta_t, double epsilon, double a, double q_3_0, double p_3_0){
  int n_step = (int)(big_T/delta_t);
  //Das tiempo
  double*t;
  t = crearArregloCero(n_step);
  //Los arreglitos
  double*q_1;
  double*p_1;
  double*q_3;
  double*p_3;
  double*energia;

  //Crea los arreglos de qs y ps
  q_1 = crearArregloCero(n_step);
  p_1 = crearArregloCero(n_step);
  q_3 = crearArregloCero(n_step);
  p_3 = crearArregloCero(n_step);
  energia = crearArregloCero(n_step);

  //Condiciones Iniciales
  //Todo - Aun falta como hacerlo. Lo mas practico probablemente sera traerlas desde un archivo.
  //Por ahora sigamos las del paper para la masa 1 (a,0) y pongamos la masa 3 en (0,0)
  q_1[0] = a;
  p_1[0] = 0.0;
  q_3[0] = q_3_0;
  p_3[0] = p_3_0;

  //Primera iteracion
  double*results_1_0;
  double*results_3_0;
  results_1_0 = RK4_step_1(delta_t, q_1[0], p_1[0], q_3[0], p_3[0], epsilon);
  results_3_0 = RK4_step_3(delta_t, q_1[0], p_1[0], q_3[0], p_3[0], epsilon);
  q_1[1] = results_1_0[0];
  p_1[1] = results_1_0[1];
  q_3[1] = results_3_0[0];
  p_3[1] = results_3_0[1];
  t[1] = delta_t;
  energia[1] = hamiltoniano_t(q_1[1], p_1[1], q_3[1], p_3[1], epsilon);
  //RK4

  int i;

  for(i = 1; i<n_step;i++){
    double*results_1_0;
    double*results_3_0;
    results_1_0 = RK4_step_1(delta_t, q_1[i-1], p_1[i-1], q_3[i-1], p_3[i-1], epsilon);
    results_3_0 = RK4_step_3(delta_t, q_1[i-1], p_1[i-1], q_3[i-1], p_3[i-1], epsilon);
    q_1[i] = results_1_0[0];
    p_1[i] = results_1_0[1];
    q_3[i] = results_3_0[0];
    p_3[i] = results_3_0[1];
    t[i] = t[i-1] + delta_t;
    energia[i] = hamiltoniano_t(q_1[i], p_1[i], q_3[i], p_3[i], epsilon);
  }

  //Ya termino la integracion. Ahora imprimamos los resultados. Son 4 arreglos.
  escribirArreglos_RK4(t, q_1, q_3, p_1, p_3, energia, n_step);
}

/*
 * Das RK4 paso, el de toda la vida
 * Paso para la masa 1. 
 */
double* RK4_step_1(double delta_t, double q_1, double p_1, double q_3, double p_3, double epsilon){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double k_step, l_step;
  double *retorno;

  k1 = q_punto(p_1);
  l1 = p_punto_1(q_1, q_3, epsilon);

  k2 = q_punto( p_1 + l1*delta_t*0.5);
  l2 = p_punto_1( q_1 + k1*delta_t*0.5, q_3, epsilon);

  k3 = q_punto( p_1 + l2*delta_t*0.5);
  l3 = p_punto_1( q_1 + k2*delta_t*0.5, q_3, epsilon);

  k4 = q_punto( p_1 + l3*delta_t);
  l4 = p_punto_1( q_1 + k3*delta_t, q_3, epsilon);

  k_step = (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  l_step = (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;

  retorno = crearArregloCero(2);
  retorno[0] = q_1+k_step;
  retorno[1] = p_1+l_step;

  return retorno;
}

/*
 * Das RK4 paso, el de toda la vida
 * Paso para la masa 3. 
 */
double* RK4_step_3(double delta_t, double q_1, double p_1, double q_3, double p_3, double epsilon){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double k_step, l_step;
  double *retorno;

  k1 = q_punto(p_3);
  l1 = p_punto_3(q_1, q_3, epsilon);

  k2 = q_punto( p_3 + l1*delta_t*0.5);
  l2 = p_punto_3( q_1, q_3 + k1*delta_t*0.5, epsilon);

  k3 = q_punto( p_3 + l2*delta_t*0.5);
  l3 = p_punto_3( q_1, q_3 + k2*delta_t*0.5, epsilon);

  k4 = q_punto( p_3 + l3*delta_t);
  l4 = p_punto_3( q_1, q_3 + k3*delta_t, epsilon);

  k_step = (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  l_step = (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;

  retorno = crearArregloCero(2);
  retorno[0] = q_3+k_step;
  retorno[1] = p_3+l_step;

  return retorno;
}

/*
 * Funcion que escribe 4 arreglos TODOS de tamano n_step (dado su tamano).
 * En el paper nos dicen que solo queremos los puntos q_3 y p_3 cuando p_1 es 0 (o pasa por 0).
 * Queremos la energia para mostrar conservacion de energia en la ultima grafica
 * Asi ademas no tenemos que escribir tantas cosas.
 */
void escribirArreglos_RK4(double*t, double*q_1, double*q_3, double*p_1, double*p_3, double* energia, int n_step)
{
  FILE* archivo;
  //Sabemos que solo hay 6 arreglos
  //Sabemos que la condicion inicial es la posicion inicial de la masa grande, o sea de la masa 1, o sea q1

  char bufX[20];
  int i;
  sprintf(bufX, "%f", q_1[0]);
  char n1[50], n3[50];
   strcpy(n1,  "int_RK4_");
   strcpy(n3, ".dat");
   strcat(n1, bufX);
   strcat(n1, n3);
   archivo = fopen(n1, "a");
   //Aca verifica que valga la pena escribir el primer de punto
   //Verifica si la velocidad de la masa 1 es 0
   if(p_1[0]==0){
     fprintf(archivo, "%f \t %.15e \t %.15e \t %.15e \n", t[0], q_3[0], p_3[0], energia[0]);
   }
  for(i=1;i<n_step;i++){
    //Aca verifica que valga la pena escribir el resto de puntos
    //Verifica si la velocidad de la masa 1 es 0 o si acaba de pasar por 0 (o sea si cambio de signo)
    if(p_1[i]==0 || (p_1[i-1]<0 && p_1[i]>0) || (p_1[i]<0 && p_1[i-1]>0) ){
      fprintf(archivo, "%f \t %.15e \t %.15e \t %.15e \n", t[i], q_3[i], p_3[i], energia[i]);
    }    
  }
  //fclose(archivo);
}

/*
 *---------------------------------------------------------------------------------------------
 * Integracion Simpletica
 */

/*
 * Metodo que realiza toda la integracion simpletica
 */
void Simpletic_Integration(double big_T, double delta_t, double epsilon, double a, double q_3_0, double p_3_0){
  int n_step = (int)(big_T/delta_t);
  //Das tiempo
  double*t;
  t = crearArregloCero(n_step);
  //Los arreglitos
  double*q_1;
  double*p_1;
  double*q_3;
  double*p_3;
  double*energia;

  //Crea los arreglos de qs y ps
  q_1 = crearArregloCero(n_step);
  p_1 = crearArregloCero(n_step);
  q_3 = crearArregloCero(n_step);
  p_3 = crearArregloCero(n_step);
  energia = crearArregloCero(n_step);

  //Les alphas
  double alpha_0, alpha_1;
  alpha_0 = -(pow(2.0,(1.0/3.0)))/(2.0-pow(2,(1.0/3.0)));
  alpha_1 = 1.0/(2.0-pow(2,(1.0/3.0)));

  //Condiciones Iniciales
  //Todo - Aun falta como hacerlo. Lo mas practico probablemente sera traerlas desde un archivo.
  //Por ahora sigamos las del paper para la masa 1 (a,0) y pongamos la masa 3 en (0,0)
  q_1[0] = a;
  p_1[0] = 0.0;
  q_3[0] = q_3_0;
  p_3[0] = p_3_0;

  //Primera iteracion
  double*results_0;
  results_0 = Simpletic_Integration_step( alpha_1, alpha_0, q_1[0], p_1[0], q_3[0], p_3[0], epsilon);
  q_1[1] = results_0[0];
  p_1[1] = results_0[1];
  q_3[1] = results_0[2];
  p_3[1] = results_0[3];
  t[1] = delta_t;
  energia[1] = hamiltoniano_t(q_1[1], p_1[1], q_3[1], p_3[1], epsilon);
  //Integracion Simpletica

  int i;

  for(i = 1; i<n_step;i++){
    double* results;
    results = Simpletic_Integration_step( alpha_1, alpha_0, q_1[i-1], p_1[i-1], q_3[i-1], p_3[i-1], epsilon);
    q_1[i] = results[0];
    p_1[i] = results[1];
    q_3[i] = results[2];
    p_3[i] = results[3];
    t[i] = t[i-1]+delta_t;
    energia[i] = hamiltoniano_t(q_1[i], p_1[i], q_3[i], p_3[i], epsilon);
  }

  //Ya termino la integracion. Ahora imprimamos los resultados. Son 6 arreglos.
  escribirArreglos_simpletico(t, q_1, q_3, p_1, p_3, energia, n_step);
}

/*
 * Das big paso simpletico. 
 */
double* Simpletic_Integration_step(double alpha_1, double alpha_0, double q_1, double p_1, double q_3, double p_3, double epsilon){
  double* dasBigReturn;
  //El plan es hacer un lambda gigante con alpha_1 luego otro con alpha_0 y luego otro con alpha_1
  //Primero hacemos todo lo de la masa 1 porque p_3 y q_3  no salen en eso.
  double*return_1_1;
  double*return_1_2;
  double*return_1_3;
  return_1_1 = labmda_mayuscula_1(alpha_1, q_1, q_3, p_3, epsilon); 
  q_1 = return_1_1[0];
  p_1 = return_1_1[1];
  return_1_2 = labmda_mayuscula_1(alpha_0, q_1, q_3, p_3, epsilon);
  q_1 = return_1_2[0];
  p_1 = return_1_2[1];
  return_1_3 = labmda_mayuscula_1(alpha_1, q_1, q_3, p_3, epsilon);
  q_1 = return_1_3[0];
  p_1 = return_1_3[1];
  //Pasamos ahora a hacer lo de la masa 3
  double*return_3_1;
  double*return_3_2;
  double*return_3_3;
  return_3_1 = labmda_mayuscula_3(alpha_1, q_1, q_3, p_3, epsilon); 
  q_3 = return_3_1[0];
  p_3 = return_3_1[1];
  return_3_2 = labmda_mayuscula_3(alpha_0, q_1, q_3, p_3, epsilon);
  q_3 = return_3_2[0];
  p_3 = return_3_2[1];
  return_3_3 = labmda_mayuscula_3(alpha_1, q_1, q_3, p_3, epsilon);
  q_3 = return_3_3[0];
  p_3 = return_3_3[1];
  //Das Big Return
  dasBigReturn = crearArregloCero(4);
  dasBigReturn[0] = q_1;
  dasBigReturn[1] = p_1;
  dasBigReturn[2] = q_3;
  dasBigReturn[3] = p_3;
  return dasBigReturn;
}

/*
 * Das little paso simpletico
 * Para la masa 3
 */
double* labmda_mayuscula_3(double alpha,double q_1, double q_3, double p_3, double epsilon){
  //Lambda gigante es: gamma(u,tau/2) luego gamma(a,tau) luego gamma(u, tau/2)
  double*return_1;
  double*return_2;
  double*return_3;
  double*bigReturn;
  bigReturn = crearArregloCero(2);
  return_1 = gamma_U_3(alpha/2.0, q_1, q_3, p_3, epsilon);
  q_3 = return_1[0];
  p_3 = return_1[1];
  return_2 = gamma_A_3(alpha, q_1, q_3, p_3, epsilon);
  q_3 = return_2[0];
  p_3 = return_2[1];
  return_3 = gamma_U_3(alpha/2.0, q_1, q_3, p_3, epsilon);
  q_3 = return_3[0];
  p_3 = return_3[1];
  bigReturn[0]=q_3;
  bigReturn[1]=p_3;
  return bigReturn;
}

/*
 * Das little paso simpletico
 * Para la masa 1
 */
double* labmda_mayuscula_1(double alpha,double q_1, double q_3, double p_1, double epsilon){
  //Lambda gigante es: gamma(u,tau/2) luego gamma(a,tau) luego gamma(u, tau/2)
  double*return_1;
  double*return_2;
  double*return_3;
  double*bigReturn;
  bigReturn = crearArregloCero(2);
  return_1 = gamma_U_1(alpha/2.0, q_1, q_3, p_1, epsilon);
  q_1 = return_1[0];
  p_1 = return_1[1];
  return_2 = gamma_A_1(alpha, q_1, q_3, p_1, epsilon);
  q_1 = return_2[0];
  p_1 = return_2[1];
  return_3 = gamma_U_1(alpha/2.0, q_1, q_3, p_1, epsilon);
  q_1 = return_3[0];
  p_1 = return_3[1];
  bigReturn[0]=q_1;
  bigReturn[1]=p_1;
  return bigReturn;
}


/*
 * -----------------------------------------------------------------------------
 * Das mini paso simpletico
 * Los A_prime se refieren a q puntico y los U_prime a - p puntico
 * Los A son 1/2 p*p y los U 1/2 q*q
 * gamma_U(q,p) = (q ,p - tau*Uprime) = (q ,p + tau*p_puntico)
 * gamma_A(q,p) = (q + tau*Aprime, p) = (q + tau*q_puntico, p)
 * El alpha juega el rol del tau
 * Ya todo empieza a ser claro
 * si la opcion es 1 es un U y si la opcion es 0 es un A
 * en vez de jugar con opciones, mejor un metodo gamma por opcion.
 * ----------------------------------------------------------------------------
 */

/*
 * Das mini paso simpletico
 * Gamma U para p1
 */
double* gamma_U_1 (double alpha, double q_1, double q_3, double p, double epsilon){
  double *retorno;
  double new_q, new_p;
  retorno = crearArregloCero(2);
  new_q = q_1;
  new_p = p + alpha*p_punto_1(q_1,  q_3, epsilon);
  retorno[0]=new_q;
  retorno[1]=new_p;
  return retorno;
}

/*
 * Das mini paso simpletico
 * Gamma U para p3
 */
double* gamma_U_3 (double alpha, double q_1, double q_3, double p, double epsilon){
  double new_q, new_p;
  double *retorno;
  retorno = crearArregloCero(2);
  new_q = q_3;
  new_p = p + alpha*p_punto_3(q_1,  q_3, epsilon);
  retorno[0]=new_q;
  retorno[1]=new_p;
  return retorno;
}


/*
 * Das mini paso simpletico
 * Gamma A para q1
 */
double* gamma_A_1 (double alpha, double q_1, double q_3, double p, double epsilon){
  double new_q, new_p;
  double *retorno;
  retorno = crearArregloCero(2);
  new_q = q_1 + alpha*q_punto(p);
  new_p = p;
  retorno[0]=new_q;
  retorno[1]=new_p;
  return retorno;
}

/*
 * Das mini paso simpletico
 * Gamma A para q3
 */
double* gamma_A_3 (double alpha, double q_1, double q_3, double p, double epsilon){
  double new_q, new_p;
  double *retorno;
  retorno = crearArregloCero(2);
  new_q = q_3 + alpha*q_punto(p);
  new_p = p;
  retorno[0]=new_q;
  retorno[1]=new_p;
  return retorno;
}

/*---------------------------------------------------------------------------
 *Funciones que se usan en ambos metodos
 *---------------------------------------------------------------------------
/*

/*
 * q punto con m_3=0, m_1*q_1+m_2*q_2=0, m_1=m_2 y q_1=-q_2 y p_1=-p_2
 */
double q_punto(double p){
  return p;
}

/*
 * Devuelve el cuadrado de un numero
 */
double square(double x){
  return x*x;
}

/*
 * p_1 punto con m_3=0, m_1*q_1+m_2*q_2=0, m_1=m_2 y q_1=-q_2 y p_1=-p_2
 */
double p_punto_1(double q_1,  double q_3, double epsilon){
  double nom, denom;
  nom = -2.0*q_1;
  denom = pow((square(2.0*q_1)+square(epsilon)),1.5);
  return (nom/denom);
}

/*
 * p_3 punto con m_3=0, m_1*q_1+m_2*q_2=0, m_1=m_2 y q_1=-q_2 y p_1=-p_2
 */
double p_punto_3(double q_1,double q_3, double epsilon){
  double primer_nom, primer_denom, segundo_nom, segundo_denom;
  primer_nom = (q_1-q_3);
  primer_denom = pow((square(q_1-q_3)+square(epsilon/2.0)),1.5);
  segundo_nom = (q_1+q_3);
  segundo_denom = pow((square(q_1-q_3)+square(epsilon/2.0)),1.5);
  return (primer_nom/primer_denom) - (segundo_nom/segundo_denom);
}

/*
 *Hamiltoniano del sistema q_1,p_1 con m_3=0, m_1*q_1+m_2*q_2=0, m_1=m_2 y q_1=-q_2 y p_1=-p_2
 */
double hamiltoniano_1(double q_1, double p_1, double epsilon){
  double primer_nom, primer_denom, segundo_nom, segundo_denom;
  primer_nom = square(p_1);
  primer_denom = 2.0;
  segundo_nom = 1;
  segundo_denom = 2.0*pow(square(2.0*q_1)+square(epsilon),0.5);
  return (primer_nom/primer_denom)-(segundo_nom/segundo_denom);
}

/*
 *Hamiltoniano del sistema q_3,p_3 con m_3=0, m_1*q_1+m_2*q_2=0, m_1=m_2 y q_1=-q_2 y p_1=-p_2
 */
double hamiltoniano_3(double q_3, double p_3, double epsilon){
  double primer_nom, primer_denom, segundo_nom, segundo_denom;
  primer_nom = square(p_3);
  primer_denom = 2.0;
  segundo_nom = 1;
  segundo_denom = 2.0*pow(square(2.0*q_3)+square(epsilon),0.5);
  return (primer_nom/primer_denom)-(segundo_nom/segundo_denom);
}

/*
 *Hamiltoniano total del sistema con m_3=0, m_1*q_1+m_2*q_2=0, m_1=m_2 y q_1=-q_2 y p_1=-p_2
 */
double hamiltoniano_t(double q_1, double p_1,double q_3, double p_3, double epsilon){
  return hamiltoniano_3(q_3, p_3, epsilon) + hamiltoniano_1(q_1, p_1, epsilon); 
}


/*
 *Funcion que crea y retorna un vector de n_points posiciones con 0.0 en ellas.
 */
double*crearArregloCero(int n_points)
{
  int i;
  double* arregloRespuesta;
  arregloRespuesta= malloc(n_points * sizeof(double));
  //Inicialicemoslo en 0
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=0.0;
    }
  return arregloRespuesta;
}

/*
 * Funcion que escribe 4 arreglos TODOS de tamano n_step (dado su tamano).
 * En el paper nos dicen que solo queremos los puntos q_3 y p_3 cuando p_1 es 0 (o pasa por 0).
 * Queremos la energia para mostrar conservacion de energia en la ultima grafica
 * Asi ademas no tenemos que escribir tantas cosas.
 */
void escribirArreglos_simpletico(double*t, double*q_1, double*q_3, double*p_1, double*p_3, double* energia, int n_step)
{
  FILE* archivo;
  //Sabemos que solo hay 6 arreglos
  //Sabemos que la condicion inicial es la posicion inicial de la masa grande, o sea de la masa 1, o sea q1
  char bufX[20];
  int i;
  sprintf(bufX, "%f", q_1[0]);
  char n1[50], n3[50];
   strcpy(n1,  "int_simpletico_");
   strcpy(n3, ".dat");
   strcat(n1, bufX);
   strcat(n1, n3);
   archivo = fopen(n1, "a");
   //Aca verifica que valga la pena escribir el primer de punto
   //Verifica si la velocidad de la masa 1 es 0
   if(p_1[0]==0){
     fprintf(archivo, "%f \t %.15e \t %.15e \t %.15e \n", t[0], q_3[0], p_3[0], energia[0]);
   }
  for(i=1;i<n_step;i++){
    //Aca verifica que valga la pena escribir el resto de puntos
    //Verifica si la velocidad de la masa 1 es 0 o si acaba de pasar por 0 (o sea si cambio de signo)
    if(p_1[i]==0 || (p_1[i-1]<0 && p_1[i]>0) || (p_1[i]<0 && p_1[i-1]>0) ){
      fprintf(archivo, "%f \t %.15e \t %.15e \t %.15e \n", t[i], q_3[i], p_3[i], energia[i]);
    }    
  }
  fclose(archivo);
}


/*---------------------------------------------------------------------------
 *Funciones que no estoy usando pero pueden llegar a ser utiles
 *---------------------------------------------------------------------------
/*

/*
 * Metodo que crea una matriz n*m con 0 en todas las entradas.
 * @params n numero de filas
 * @params m numero de columnas
 * @returns *double matriz de ceros n*m
 */
double *crear_matriz( int n, int m){
  double *matrix;
  int i;
  int j;
  matrix = malloc(n * m * sizeof(double));

  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
       matrix[i*m + j]=0.0;
    }
  }    
  return matrix;
}

/*
 *Funcion que crea y retorna un vector de n_points puntos equiespaciados entre x_ini y x_fin
 */
double*crearArregloEquiEspaciado(double x_ini,double x_fin,int n_points)
{
  int i;
  double* arregloRespuesta;
  arregloRespuesta= malloc(n_points * sizeof(double));
  double paso;
  paso=(x_fin-x_ini)/n_points;
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=x_ini+paso*i;
    }
  return arregloRespuesta;
}


/*
 *Funcion que dado un arreglo y su tamano, crea y retorna una copia de este.
 */
double*copiarArreglo(double*original, int n_points)
{
  int i;
  double* arregloRespuesta;
  arregloRespuesta=malloc(n_points * sizeof(double));
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=original[i];
    }
  return arregloRespuesta;
}
/*
 *Funcion que dada una matriz y su tamano, crea y retorna una copia de esta.
 */
double*copiarMatriz(double*original, int n, int m)
{
double *matricitaResultado=crear_matriz( n, m);
 int i,j;
  for (i=0;i<n;i++)
    {
    for(j=0;j<m;j++)
      {
	matricitaResultado[i*m+j]=0.0;//Ponemos la matriz en ceros para evitar problemas
	matricitaResultado[i*m+j]=original[i*m+j];
      }
    }
  return matricitaResultado;
}
