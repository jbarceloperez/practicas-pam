// mulmatcua.c
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h> 

int THREADS = 4;
#define TAM 16
#define GRANO 1

///////////////////////////////////////////////////////////////

void mv(double *m1, double *m2, double *mr,int t,int chunk) {
  int i,j,k,iam,nprocs;
  double suma;
  #pragma omp parallel private(iam,nprocs) 
  {
    nprocs=omp_get_num_threads();
    iam=omp_get_thread_num();
    #pragma omp master
      THREADS=nprocs;
    #pragma omp for private(i,suma,j) schedule(dynamic,chunk)
    for (i = 0; i < t; i++) {
        #ifdef DEBUG
            printf("thread %d fila %d \n",iam,i);
        #endif
        suma=0.;
        for(j=0;j<t;j++){
            for(k=0;k<t;k++)
                suma += m1[i*t+k] * m2[j+k*t];
            mr[i*t+j]=suma;
            suma = 0.;
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
void escribir(double *m, int t) {
  int i,j;
  for (i=0;i<t;i++) {
    for(j=0;j<t;j++)
      printf("%.4lf ",m[i*t+j]);
    printf("\n");
  }
}
///////////////////////////////////////////////////////////////

int main(int argc,char *argv[]) {

    double *m1, *m2, *mr;   // matrices 1 y 2 a multiplicar, y mr es la matriz donde se guardará el resultado
    int opt;
    int t = TAM;        // tam de la matriz cuadrada
    int chunk = GRANO;  // cuantas filas consecutivas le tocará a cada thread, para el schedule()
    double start,fin,tiempo,Mflops;

    optind = 1;
    while ((opt = getopt(argc, argv, "t:c:h")) != -1)
    {
        switch (opt)
        {
        case 't':
            t = atoi(optarg);
            break;
        
        case 'c':
            chunk = atoi(optarg);
            break;
        
        case 'h':
            printf("\n\n USO %s [-t TAM_MATRIZ] [-c TAM_SCHEDULE] \n\n",argv[0]);
            return -1;
        }
    }


    // inicializar matrices
    m1 = (double *) malloc(sizeof(double)*t*t); //matriz 1
    m2 = (double *) malloc(sizeof(double)*t*t); //matriz 2
    mr = (double *) malloc(sizeof(double)*t*t); //matriz resultado
    initializealea(m1,t*t);
    initializealea(m2,t*t);

   
    printf("m1:\n");
    escribir(m1,t);
    printf("m2:\n");
    escribir(m2,t);
 

    start=omp_get_wtime();
    mv(m1, m2, mr, t, chunk);
    fin=omp_get_wtime();
    tiempo=fin-start;
    if(tiempo==0.)
    {
      printf("No hay suficiente precision\n");
    }
    else 
    {
      Mflops=((2.*t*t)/tiempo)/1000000.;
      printf("\n  Threads %d, tamano %d\n    segundos: %.6lf, Mflops: %.6lf, Mflops por thread: %.6lf\n",THREADS,t,tiempo,Mflops,Mflops/THREADS);
    }


    printf("m resultado:\n");
    escribir(mr,t);//debug

    free(m1);
    free(m2);
    free(mr);
}
