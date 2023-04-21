#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h> 


///////////////////////////////////////////////////////////////
void mmcua(double *a, double *b, double *c, int n, int f) {
  int i,j,k,iam,nprocs;
  double suma;

  #pragma omp parallel private(iam,nprocs)
  {
    nprocs=omp_get_num_threads();
    iam=omp_get_thread_num();

    #pragma omp for private(i,suma,j) schedule(dynamic,f)
    for (i=0; i<n; i++) {
      
      #ifdef DEBUG
        printf("thread %d fila %d \n",iam,i);
      #endif
      
      for(j=0; j<n; j++) {
        suma = 0.;
        for (k=0; k<n; k++){
          suma += a[i*n+k] * b[j*n+k];
        }
        c[i*n+j] = suma;
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
  int n,f,t;
  double start,fin,tiempo,Mflops;
  double *a,*b,*c;
  
  if (argc!=4)
  {
    printf("\n\n USO %s <dim_mat_n> <num_filas_consecutivas_F> <num_hil_t> \n\n",argv[0]);
    return -1;
  }
  
  n=atoi(argv[1]);
  f=atoi(argv[2]);
  t=atoi(argv[3]);
  omp_set_num_threads(t);

  // inicializar matrices
  a = (double *) malloc(sizeof(double)*n*n);
  b = (double *) malloc(sizeof(double)*n*n);
  c = (double *) malloc(sizeof(double)*n*n); //matriz resultado
  initializealea(a,n*n);
  initializealea(b,n*n);

  #ifdef DEBUG
    printf("a:\n");
    escribir(a,n,n,n);
    printf("b:\n");
    escribir(b,n,n,n);
  #endif
  
  start=omp_get_wtime();
  mmcua(a, b, c, n, f);
  fin=omp_get_wtime();
  tiempo=fin-start;
  if(tiempo==0.)
  {
    printf("No hay suficiente precision\n");
  }
  else 
  {
    Mflops=((2.*n*n)/tiempo)/1000000.;
    printf("\n  Threads %d, tamano %d\n    segundos: %.6lf, Mflops: %.6lf, Mflops por thread: %.6lf\n",t,n,tiempo,Mflops,Mflops/t);
  }

  #ifdef DEBUG
    printf("resultado:\n");
    escribir(c,n,n,n);//debug
  #endif

  free(a);
  free(b);
  free(c);

}
