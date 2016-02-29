#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Definicion de constantes
// Opciones para ejecutar el codigo
// La opcion depende del valor de a
#define FIG_9_10_11 1
#define FIG_12_13_14 2
#define FIG_15 3
#define FIG_16 4

double* crearArregloCero(int n_points);
double* crearArregloEquiEspaciado(double x_ini,double x_fin,int n_points);
void escribirArreglos(double*q_3, double*p_3, int n_condiciones, int opcion);


/*
 * Metodo que se ejecuta con el programa
 * Se espera que se ejecute con ./Das_ICs.x OPTION_EXECUTOR(1, 2, 3 o 4)
 */
int main(int argc, char **argv){
  int opcion;
  double *q_3;
  double *p_3;
  opcion = atoi(argv[1]);

  if(opcion == FIG_9_10_11){
    //Para las fig 9: Los limites de la figura que queremos son: -4 a 4 para p_3 y q_3
    //Con alto detalle.
    int n_condiciones=120;
    
    double q_3_max = 4.0;
    double q_3_min = -4.0;
    double p_3_max = 3.0;
    double p_3_min = -3.0;
    q_3 = crearArregloEquiEspaciado(q_3_min, q_3_max, n_condiciones);
    p_3 = crearArregloEquiEspaciado(p_3_min, p_3_max, n_condiciones);

    escribirArreglos(q_3, p_3, n_condiciones, opcion);
   }
  else if (opcion == FIG_12_13_14){
    //Para las figuras 12_13_14 alrededor de: -2.5 a 2.5 para p_3 y -2 a 2 para q_3 con mas detalle alrededor de: p_3 [-1.5;1.5] q_3 [-0.9990;0.990] y p_3 [-0.1;0.1] q_3 [-0.01;0.01]
    int n_condiciones=120;
    
    double q_3_max = 2.5;
    double q_3_min = -2.5;
    double p_3_max = 3.0;
    double p_3_min = -3.0;
    q_3 = crearArregloEquiEspaciado(q_3_min, q_3_max, n_condiciones);
    p_3 = crearArregloEquiEspaciado(p_3_min, p_3_max, n_condiciones);

    escribirArreglos(q_3, p_3, n_condiciones, opcion);
  }
  else if(opcion == FIG_15){
    //Para las figuras 15 alrededor de: p_3 [-2.5;2.5] q_3 [-3;3]
    int n_condiciones=120;
    
    double q_3_max = 3.5;
    double q_3_min = -3.5;
    double p_3_max = 3.0;
    double p_3_min = -3.0;
    q_3 = crearArregloEquiEspaciado(q_3_min, q_3_max, n_condiciones);
    p_3 = crearArregloEquiEspaciado(p_3_min, p_3_max, n_condiciones);

    escribirArreglos(q_3, p_3, n_condiciones, opcion);
  }
  else if(opcion == FIG_16){
    //Para las figuras 16 alrededor de: p_3 [-2.5;2.5] q_3 [-3;3]
    int n_condiciones=120;
    
    double q_3_max = 3.5;
    double q_3_min = -3.5;
    double p_3_max = 3.0;
    double p_3_min = -3.0;
    q_3 = crearArregloEquiEspaciado(q_3_min, q_3_max, n_condiciones);
    p_3 = crearArregloEquiEspaciado(p_3_min, p_3_max, n_condiciones);

    escribirArreglos(q_3, p_3, n_condiciones, opcion);
  }

  return 0;
}

void escribirArreglos(double*q_3, double*p_3, int n_condiciones, int opcion)
{
  FILE* archivo;
  //Sabemos que solo hay 6 arreglos
  //Sabemos que la condicion inicial es la posicion inicial de la masa grande, o sea de la masa 1, o sea q1

  char bufX[20];
  int i;
  sprintf(bufX, "%d", opcion);
  char n1[50], n3[50];
  if(opcion == FIG_9_10_11){
   strcpy(n1,  "ics_fig_9_10_11");
  }
  else if(opcion == FIG_12_13_14){
   strcpy(n1,  "ics_fig_12_13_14");
  }
  else if(opcion == FIG_15){
   strcpy(n1,  "ics_fig_15");
  }
  else if(opcion == FIG_16){
   strcpy(n1,  "ics_16");
  }
   strcpy(n3, ".dat");
   strcat(n1, bufX);
   strcat(n1, n3);
   archivo = fopen(n1, "w");

   fprintf(archivo, "%d  \n", n_condiciones);

  for(i=0;i<n_condiciones;i++){
      fprintf(archivo, " %.15e \t %.15e\n", q_3[i], p_3[i]);
  }    
  fclose(archivo);
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
