#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h> 


/**
 * Calculo del producto de dos matrices AxB, siendo A m x km y B km x n, mediante
 * un esquema por bloques (tiling) de tamaño b, y almacena el resultado en C.
 * @param a matriz A
 * @param b matriz B, indexada por columnas
 * @param c matriz C
 * @param m num de filas de A y C 
 * @param n num de columnas de B y C
 * @param km num de filas de B y de columnas de A
 * @param bloque tamano del bloque
*/
void seqmmblo(double *a, double *b, double *c, int m, int n, int km, int bloque) {
    int i,j,k,iout,jout,kout;

    // estos dos bucles (iout y jout) for recorren los diferentes bloques de la
    // matriz resultado que se deben calcular
    for (iout=0; iout<m; iout+=bloque) {      
        for(jout=0; jout<n; jout+=bloque) {
            // luego el bucle de kout calcula para cada bloque de la matriz resultado
            // los bloques de a y b que debe recorrer y operar
            for (kout=0; kout<km; kout+=bloque) {
                // CALCULO DE UN BLOQUE DE C, se usa el mismo algoritmo de seqmulmatcuablo.c
                // pero usando km para las columnas de A y filas de B, para indexar los elementos
                // de la matriz resultado en el calculo
                for (i=iout; i<iout+bloque; i++){
                    for (j=jout; j<jout+bloque; j++) {
                        for (k=kout; k<kout+bloque; k++){
                            c[i*n+j] += a[i*km+k] * b[j*km+k];    // Matriz B indexada por columnas, mejorar localidad espacial
                        }
                    }
                }
                // FIN CALCULO DE UN BLOQUE DE C
            }    
        }
    }
}

///////////////////////////////////////////////////////////////
void initialize(double *m,int t) {
  int i;
  for (i = 0; i < t; i++)
    m[i] = (double)(i);
}

///////////////////////////////////////////////////////////////
void initializealea(double *m, int t) {
  int i;
  for (i = 0; i < t; i++)
    m[i] = (double)rand()/RAND_MAX;
}
///////////////////////////////////////////////////////////////
void escribir(double *m, int fm,int cm) {
  int i,j;
  for (i = 0; i < fm; i++) {
    for(j=0;j<cm;j++)
      printf("%.4lf ",m[i*cm+j]);
    printf("\n");
  }
}
///////////////////////////////////////////////////////////////
void escribir_col(double *m, int fm,int cm) {
  int i,j;
  for (i = 0; i < fm; i++) {
    for(j=0;j<cm;j++)
      printf("%.4lf ",m[i+fm*j]);
    printf("\n");
  }
}

///////////////////////////////////////////////////////////////
int main(int argc,char *argv[]) {
  int m,n,k,b;
  double start,fin,tiempo,Mflops;
  double *A,*B,*C;
  
  if (argc!=5)
  {
    printf("\n\n USO %s <dim_mat_m> <dim_mat_n> <dim_mat_k> <tam_blo_b> \n\n",argv[0]);
    return -1;
  }
  
  m=atoi(argv[1]);
  n=atoi(argv[2]);
  k=atoi(argv[3]);
  b=atoi(argv[4]);

  // inicializar matrices
  A = (double *) malloc(sizeof(double)*m*k);
  B = (double *) malloc(sizeof(double)*k*n);    
  C = (double *) calloc(m*n, sizeof(double));   //matriz resultado, inicializada a cero
  initializealea(A,m*k);
  initializealea(B,k*n);    // la matriz b se interpreta como si estuviera almacenada por columnas, como son valores aleatorios da igual, en otro caso haría falta transponerla

  #ifdef DEBUG
    printf("a:\n");
    escribir(A,m,k);
    printf("b:\n");
    escribir_col(B,k,n);
  #endif
  
  start=omp_get_wtime();
  seqmmblo(A, B, C, m, n, k, b);
  fin=omp_get_wtime();
  tiempo=fin-start;
  if(tiempo==0.)
  {
    printf("No hay suficiente precision\n");
  }
  else 
  {
    Mflops=((2.*n*n)/tiempo)/1000000.;
    printf("\n  m=%d, n=%d, k=%d, tam blo: %d\n    segundos: %.6lf, Mflops: %.6lf, Mflops por thread: %.6lf\n",m,n,k,b,tiempo,Mflops,Mflops/1);
  }

  #ifdef DEBUG
    printf("resultado:\n");
    escribir(C,m,n);//debug
  #endif

  free(A);
  free(B);
  free(C);

}
