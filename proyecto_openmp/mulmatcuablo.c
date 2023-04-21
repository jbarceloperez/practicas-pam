#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h> 

/**
 * Calculo del producto de dos matrices cuadradas AxB, ambas nxn, mediante
 * un esquema por bloques (tiling) de tamaño b, y almacena el resultado en C.
 * Usa directrices OpenMP para paralelizar los bucles for.
 * @param a matriz A
 * @param b matriz B, indexada por columnas
 * @param c matriz C
 * @param n tamano de la matriz
 * @param bloque tamano del bloque
*/
void mmcuablo(double *a, double *b, double *c, int n, int bloque) {
    int i,j,k,iout,jout,kout,iam;

    #pragma omp parallel private(iam)
    {
        iam = omp_get_thread_num();

        // estos dos bucles (iout y jout) for recorren los diferentes bloques de la
        // matriz resultado que se deben calcular, por eso se paraleliza con collapse(2),
        // para que OpenMP reparta cada bloque de C a un hilo distinto.
        #pragma omp for private (iout,jout,kout,i,j,k) schedule(dynamic) collapse(2)
        for (iout=0; iout<n; iout+=bloque) {      
            for(jout=0; jout<n; jout+=bloque) {
                // luego el bucle de kout calcula para cada bloque de la matriz resultado
                // los bloques de a y b que debe recorrer y operar. cada iteración la ejecuta
                // solamente un hilo.
                #ifdef DEBUG
                    printf("thread %d calculando C(%d,%d) <-- threads distintos\n",iam,iout/bloque,jout/bloque);
                #endif
                for (kout=0; kout<n; kout+=bloque) {
                    // CALCULO DE UN BLOQUE DE C
                    #ifdef DEBUG
                        printf("\tthread %d calculando C(%d,%d) <-- mismo thread\n",iam,iout/bloque,jout/bloque);
                    #endif
                    for (i=iout; i<iout+bloque; i++) {
                        for (j=jout; j<jout+bloque; j++) {
                            for (k=kout; k<kout+bloque; k++) {
                                c[i*n+j] += a[i*n+k] * b[j*n+k];    // Matriz B indexada por columnas, mejorar localidad espacial
                            }
                        }
                    } // FIN CALCULO DE UN BLOQUE DE C
                }
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
void escribir(double *m, int fm,int cm,int ldm) {
  int i,j;
  for (i = 0; i < fm; i++) {
    for(j=0;j<cm;j++)
      printf("%.4lf ",m[i*ldm+j]);
    printf("\n");
  }
}

///////////////////////////////////////////////////////////////
int main(int argc,char *argv[]) {
  int n,b,t;
  double start,fin,tiempo,Mflops;
  double *A,*B,*C;
  
  if (argc!=4)
  {
    printf("\n\n USO %s <dim_mat_n> <tam_blo_b> <num_hil_t> \n\n",argv[0]);
    return -1;
  }
  
  n=atoi(argv[1]);
  b=atoi(argv[2]);
  t=atoi(argv[3]);

  omp_set_num_threads(t);

  // inicializar matrices
  A = (double *) malloc(sizeof(double)*n*n);
  B = (double *) malloc(sizeof(double)*n*n);    
  C = (double *) calloc(n*n, sizeof(double));   //matriz resultado, inicializada a cero
  initializealea(A,n*n);
  initializealea(B,n*n);    // la matriz b se interpreta como si estuviera almacenada por columnas, como son valores aleatorios da igual, en otro caso haría falta transponerla

  #ifdef DEBUG
    printf("a:\n");
    escribir(A,n,n,n);
    printf("b:\n");
    escribir(B,n,n,n);
  #endif
  
  start=omp_get_wtime();
  mmcuablo(A, B, C, n, b);
  fin=omp_get_wtime();
  tiempo=fin-start;
  if(tiempo==0.)
  {
    printf("No hay suficiente precision\n");
  }
  else 
  {
    Mflops=((2.*n*n)/tiempo)/1000000.;
    printf("\n  Tamano: %d, tam blo: %d\n    segundos: %.6lf, Mflops: %.6lf, Mflops por thread: %.6lf\n",n,b,tiempo,Mflops,Mflops/t);
  }

  #ifdef DEBUG
    printf("resultado:\n");
    escribir(C,n,n,n);//debug
  #endif

  // Liberar memoria
  free(A);
  free(B);
  free(C);

}
