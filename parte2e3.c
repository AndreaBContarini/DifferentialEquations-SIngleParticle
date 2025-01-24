//LFC1 - Prima prova in itinere - Parte II e III
//Andrea Belli Contarini 1916927

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//Struct
typedef struct {
  int n_passi;
  double m, x0, v0, dt, T_max, g;
  double * x;
} param;

//Definizione funzioni esterne
param leggi_input (void);
double b (double x, double v, param pd);
void RK2_ATTRITO (param pd);
double stima_v0min (param pd);
double stima_v0max (param pd);
void regione (param pd);

/*-----------------------------INIZIO MAIN-----------------------------*/
int main(int argc, char *argv[]){
  param pd;
  double vmin, vmax;
  int i, Parte;
  
  //1° check: corretta scrittura su terminale
  if (argc!=6) {
    fprintf(stderr, "CORRECT USAGE: %s gamma v0 dt T_max Parte <input.dat\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  //Assegnazione parametri iniziali
  pd = leggi_input();
  pd.g =  atof(argv[1]); //g in [0,0.2]
  pd.v0 = - atof(argv[2]);
  pd.dt = atof(argv[3]);
  pd.T_max = atof(argv[4]);
  Parte = atoi(argv[5]);
  pd.n_passi = ceil(pd.T_max/pd.dt);
  
  pd.x = (double *) calloc (pd.n_passi+1, sizeof(double));
  //controllo corretta allocazione dinamica memoria
  if (pd.x == NULL) {
    fprintf(stderr, "Errore nell'allocazione\n");
    exit(EXIT_FAILURE);
  }
  
  //2° check: dt e T_max positivi e diversi da 0
  if (pd.dt<=0 || pd.T_max<=0) {
    fprintf(stderr,"ERRORE, inserire dt>0 e T_max>0\n");
    exit(EXIT_FAILURE);
  }
  
  if (Parte==2) {
    vmin = stima_v0min(pd);
    vmax = stima_v0max(pd);
    fprintf(stdout, "Velocità minima e massima (in modulo) per cui x ----> 0:\n|vmin| = %g\n|vmax| = %g\n", fabs(vmin), fabs(vmax));
  } else if (Parte==3) {
    regione(pd);
  }  else {
    fprintf(stderr, "Inserirre solo 2 o  3 per 'Parte'\n");
    exit(EXIT_FAILURE);
  }  
  free(pd.x);
  return(EXIT_SUCCESS);
}/*-----------------------------FINE MAIN-----------------------------*/

//Funzioni esterne
param leggi_input (void){
  param res;
  FILE *fp;
  double inutile;
  
  //Lettura dati da file, comando terminale: <input.dat
  fp = fopen("input.dat", "r");
  if ((fp=fopen("input.dat", "r")) == NULL) {
    fprintf(stderr, "Errore nell'apertura del file %s\n", "input.dat");
    exit(EXIT_FAILURE);
  } else {
    fscanf(fp, "%lf", &res.m);
    fscanf(fp, "%lf", &res.x0);
    fscanf(fp, "%lf", &inutile);
    fclose(fp);
  }
  return res;
}

double b (double x, double v, param pd){
  return (exp(-x/2) * log(x) * (log(x) - (4/x))) - pd.g * v;
}

void RK2_ATTRITO (param pd) {
  int i;
  double x=pd.x0, v=pd.v0, b_n, dx, dv, xold = x;
  
  pd.x[0]=pd.x0;
  for (i=1; i<=pd.n_passi; i++) {
    
    //algoritmo runge-kutta 2
    dx = v * pd.dt;
    dv = b(x,v,pd) * pd.dt;
    
    x += (v + 0.5 * dv) * pd.dt;
    
    //SALVO LA TRAIETTORIA
    pd.x[i] = x; 
    
    v += b(xold + 0.5 * dx, v, pd) * pd.dt;
    xold = x;
  }
}

double stima_v0min (param pd){
  double vmin;
  
  pd.v0 = -0.01;
  do {
    pd.v0 -= 0.01;
    RK2_ATTRITO(pd);
  } while (pd.x[pd.n_passi]>pd.x[0]);
  
  vmin = pd.v0;
  return vmin;
}

double stima_v0max (param pd){
  double vmax;
  pd.v0 = -4.00; // v varia in [0.0, 4.0], con Delta_v = 0.1
  
  do {
    pd.v0 += 0.01;
    RK2_ATTRITO(pd);
  } while (pd.x[pd.n_passi]>pd.x[0]);
  vmax = pd.v0;
  return vmax;
}

void regione (param pd){
  int i;
  double vmax, vmin;
  pd.g = 0.01;
  
  for (i=1; i<=20; i++) {
    // gamma varia in [0.00,0.20], con Delta_g = 0.01; per cui: i da 1 a 20                               
    vmin = stima_v0min(pd);
    vmax = stima_v0max(pd);
    
    fprintf(stdout, "%g %g\n%g %g\n\n", pd.g, fabs(vmin), pd.g, fabs(vmax));
    
    pd.g += 0.01;
  }
}
