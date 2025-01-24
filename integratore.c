//LFC1 - Prima prova in itinere
//Andrea Belli Contarini 1916927

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//Struct
typedef struct {
  int n_passi;
  double m, x0, v0, dt, T_max, g;
} param;

//Definizione funzioni esterne
param leggi_input (void);
double a (double x, param pd);
double check (double x, double v);
double b (double x, double v, param pd);
void VERLET_VEL (param pd);
void RK2_ATTRITO (param pd);

/*-----------------------------INIZIO MAIN-----------------------------*/
int main(int argc, char *argv[]){
  param pd;
  int i, Parte;
  
  //1° check: corretta scrittura su terminale
  if (argc!=5) {
    fprintf(stderr, "CORRECT USAGE: %s Parte v0 dt T_max <input.dat\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  //Assegnazione parametri iniziali
  pd = leggi_input();
  Parte =  atoi(argv[1]);
  pd.v0 = - atof(argv[2]);
  pd.dt = atof(argv[3]);
  pd.T_max = atof(argv[4]);
  pd.n_passi = ceil(pd.T_max/pd.dt);
  
  //2° check: dt e T_max positivi e diversi da 0
  if (pd.dt<=0 || pd.T_max<=0) {
    fprintf(stderr,"ERRORE, inserire dt>0 e T_max>0\n");
    exit(EXIT_FAILURE);
  }
  
  if (Parte==1) {
    VERLET_VEL(pd);
  } else if (Parte==2) {
    RK2_ATTRITO (pd);
  } else {
    fprintf (stderr, "Inserire solo 1 oppure 2 per la 'Parte' da studiare\n");
    exit(EXIT_FAILURE);
  }
  
  /* La parte commentata che segue è stata utilizzata una sola volta per verificare che l'algoritmo fosse effettivamente del II ordine */
  /* for (i=0; i<7; i++) {
     VERLET_VEL(pd);
     pd.n_passi *= 2;
     pd.dt /= 2;
     }*/
  return(EXIT_SUCCESS);
}/*-----------------------------FINE MAIN-----------------------------*/

//Funzioni esterne
param leggi_input (void){
  param res;
  FILE *fp;
  
  //Lettura dati da file, comando terminale: <input.dat
  fp = fopen("input.dat", "r");
  if ((fp=fopen("input.dat", "r")) == NULL) {
    fprintf(stderr, "Errore nell'apertura del file %s\n", "input.dat");
    exit(EXIT_FAILURE);
  } else {
    fscanf(fp, "%lf", &res.m);
    fscanf(fp, "%lf", &res.x0);
    fscanf(fp, "%lf", &res.g); //fattore gamma (attrito)
    fclose(fp);
  }
  return res;
}

double a (double x, param pd){
  return exp(-x/2) * log(x) * (log(x) - (4/x));
}

double check (double x, double v){
  return 0.5 * (v*v) + 2 * (log(x)*log(x)) * exp(-x/2);
}

double b (double x, double v, param pd){
  return (exp(-x/2) * log(x) * (log(x) - (4/x))) - pd.g * v;
}

void VERLET_VEL (param pd){
  int i;
  double x=pd.x0, v=pd.v0, a_n, cost0 = check(x,v), cost, D_cost;
  
  fprintf(stdout, "#INTEGRAZIONE CON ALGORITMO VERLET VELOCITA' \n#m=%.2lf, T_max=%.4lf, dt=%.4lf \n# 1:t, 2:x, 3:v, 4:D_cost, 5:dt \n0 %.14lf %.14lf\n", pd.m, pd.T_max, pd.dt, x, v);
  
  for (i=1; i<=pd.n_passi; i++) {
    
    //algoritmo Verlet velocita'
    a_n = a(x,pd);
    x += v * pd.dt + 0.5 * a_n * (pd.dt*pd.dt);
    v += 0.5 * (a_n + a(x,pd)) * pd.dt;
    
    //calcolo della costante
    cost = check(x,v);
    D_cost = (cost - cost0) / cost0;
    
    if (i!=pd.n_passi){
      fprintf(stdout, "%g %.14lf %.14lf %g\n", i*pd.dt, x, v, D_cost);
    } else {
      fprintf(stdout, "%g %.14lf %.14lf %g %g LAST\n\n", i*pd.dt, x, v, D_cost, pd.dt); //per poter usare comando |grep da terminale
    }
  }
}

void RK2_ATTRITO (param pd) {
  int i;
  double x=pd.x0, v=pd.v0, b_n, dx, dv, xold = x;
  
  fprintf(stdout, "#INTEGRAZIONE CON ALGORITMO RUNGE-KUTTA 2 \n#m=%.2lf, T_max=%.4lf, dt=%.4lf \n# 1:t, 2:x, 3:v, 4:D_cost, 5:dt \n0 %.14lf %.14lf\n", pd.m, pd.T_max, pd.dt, x, v);
  
  for (i=1; i<=pd.n_passi; i++) {
    
    //algoritmo runge-kutta 2
    dx = v * pd.dt;
    dv = b(x,v,pd) * pd.dt;
    
    x += (v + 0.5 * dv) * pd.dt;
    v += b(xold + 0.5 * dx, v, pd) * pd.dt;
    
    if (i!=pd.n_passi){
      fprintf(stdout, "%g %.14lf %.14lf\n", i*pd.dt, x, v);
    } else {
      fprintf(stdout, "%g %.14lf %.14lf %g LAST\n\n", i*pd.dt, x, v, pd.dt);
    }
    xold = x;
  }
}
