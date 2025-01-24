//LFC1 - Prima prova in itinere - Punto C e D, Parte 1
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
double a (double x, param pd);
void VERLET_VEL (param pd);
double Calcolo_Vmin(param pd, double v);
void calcola_T (param pd);
void VERLET_VEL2 (param pd, int n_passi, double dtt);

/*-----------------------------INIZIO MAIN-----------------------------*/
int main(int argc, char *argv[]){
  param pd;
  int esco, i, Parte;
  double v=-0.800, t_star;
  
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
  
  if (Parte==1) {
    /* Finché non ottengo il primo valore di v0 affinchè x=1 continuo a girare nel while diminuendo v di 0.001 (3 cifre significative) */
    while (esco!=1) {
      v-=0.001;
      esco = Calcolo_Vmin(pd, v);
      free(pd.x);
      pd.x = (double *) calloc (pd.n_passi+1, sizeof(double));
      //controllo corretta allocazione dinamica memoria
      if (pd.x == NULL) {
	fprintf(stderr, "Errore nell'allocazione\n");
	exit(EXIT_FAILURE);
      }
    }
  } else if (Parte==2) {
    // for (i=0; i<10; i++) { 
    VERLET_VEL(pd);
    calcola_T (pd);
    /*  pd.dt /= 2;
	pd.n_passi *= 2;
	} */
  } else { //uscita se errata assegnazione di 'Parte'
    fprintf (stderr, "Inserire solo 1 oppure 2 ('Parte' da studiare; 1=c, 2=d)\n");
    exit(EXIT_FAILURE);
  }
  free(pd.x);
  return(EXIT_SUCCESS);
}/*-----------------------------FINE MAIN-----------------------------*/

//Funzioni esterne
param leggi_input (void) {
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

double a (double x, param pd) {
  return exp(-x/2) * log(x) * (log(x) - (4/x));
}

void VERLET_VEL (param pd){
  int i;
  double x=pd.x0, v=pd.v0, a_n;
  
  pd.x[0]=pd.x0;
  for (i=1; i<=pd.n_passi; i++) {
    
    //algoritmo Verlet velocita'
    a_n = a(x,pd);
    x += v * pd.dt + 0.5 * a_n * (pd.dt*pd.dt);
    //SALVO LA TRAIETTORIA
    pd.x[i] = x;
    
    v += 0.5 * (a_n + a(x,pd)) * pd.dt;
  }
}

double Calcolo_Vmin(param pd, double v){
  int i, esco;
  double vmin=v, x=pd.x0, a_n;
  pd.x[0]=pd.x0;
  
  for (i=1; i<=pd.n_passi; i++) {
    //algoritmo Verlet velocita'
    a_n = a(x,pd);
    x += v * pd.dt + 0.5 * a_n * (pd.dt*pd.dt);
    //SALVO LA TRAIETTORIA
    pd.x[i] = x;
    v += 0.5 * (a_n + a(x,pd)) * pd.dt;
    if ((pd.x[i]-1) * (pd.x[i-1]-1) <= 0) {
      //stampo su output
      fprintf(stdout, "Valore di v0_min per cui x=1:\nv0_min = %g\n", vmin);
      esco=1;
    }
  }
  return esco;
}

void calcola_T (param pd) {
  int i;
  double t_star1 = -11, t_star2;
  
  /* Alla luce della traiettoria ottenibile con gnuplut per v0=1, si fa partire il ciclo for da i=2 per evitare di ottenere valori di t_star errati (ad esempio prossimi allo 0,00... se il ciclo partisse da i=1) */
  
  for (i=2; i<pd.n_passi; i++) {
    
    //Traslazione della traiettoria di -1
    
    if ((pd.x[i]-1) * (pd.x[i+1]-1) <= 0) {
      t_star2 = ((i * pd.dt * (pd.x[i+1]-1)) - ((i+1)*pd.dt * (pd.x[i]-1)))/((pd.x[i+1]-1) - (pd.x[i]-1));
      if (t_star1 != -11) {
	fprintf(stdout, "Valori di t con v0=1 per cui x=1 (dt = %g): \n\n%g %g\n\n", pd.dt, t_star1, t_star2);
      }
      t_star1=t_star2;
      i++;
    }
  }
}
