// MulMatCuaBlo.c
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h> 

#define DEFAULT_TAM 16
#define DEFAULT_BLO 4

// imprime la matriz, para debuggear
void escribir(double *m, int t) {
  int i,j;
  for (i=0;i<t;i++) {
    for(j=0;j<t;j++)
      printf("%.4lf\t",m[i*t+j]);
    printf("\n");
  }
  printf("\t__________________________\n");
}

// IMPRIME UNA SUBMATRIZ O BLOQUE, PARA DEBUGGEAR
void escribir_submatriz(double *mr, int n, int b) {
    
    int i,j;
    for (i = 0; i < b; i++)
    {
        for (j = 0; j < b; j++)
        {
            printf("%.4lf\t",mr[i*n+j]);
        }
        printf("\n");
    }
    printf("\t________________________\n");
}

// imprime el tiempo de calculo de la matriz
void imprimir_resultados(double tiempo, int b, int n, int nthreads) {
    double mflops; 

    if (tiempo==0.)
    {
        printf("No hay suficiente precision\n");
    }
    else 
    {
        mflops=((2.*n*n)/tiempo)/1000000.;
        printf("\n  Numero de hilos: %d\n  Tam de matriz: %d\n  Tam de bloque: %d\n    segundos: %.6lf, Mflops: %.6lf\n",nthreads,n,b,tiempo,mflops);
        // printf("%.6lf\n,",mflops);
    }
}

//pone la matriz entera a 0
void reset(double *m, int b) {
    int i;
    for (i = 0; i < b*b; i++)
        m[i] = 0.;
}

// inicializa un BLOQUE con valores aleatorios 
void initializealea(double *m, int t) {     
    int i;
    for (i = 0; i < t; i++)
        m[i] = (double)rand()/RAND_MAX;
}

// multiplica dos submatrices 
void mm(double *m1, double *m2, double *mr, int t, int n) {
    int i,j,k;
    double suma;
    
    #pragma omp parallel for private(i,j,suma) schedule(dynamic, 1)
    for (i = 0; i < t; i++) {
        for(j=0;j<t;j++){
            suma=mr[i*t+j]; // se obtiene el valor guardado en aux para seguir sumándoselo
            for(k=0;k<t;k++)
                suma += m1[i*n+k] * m2[j+k*n];  // se hacen saltos de n en vez de saltos de t 
            mr[i*t+j]=suma; // se guarda el valor en la posición correspondiente de la matriz auxiliar
        }
    }
}

// INSERTA LA SUBMATRIZ AUXILIAR USADA EN EL CÁLCULO DE LA MULTIPLICACION Y SUMA EN LA MATRIZ RESULTADO
void insertar_submatriz(double *mr, double *aux, int n, int b) {
    
    int i,j;
    for (i = 0; i < b; i++)
    {
        for (j = 0; j < b; j++)
        {
            mr[i*n+j] = aux[i*b+j];
        }
    }
}

void mulmatcuablo(double *m1, double *m2, double *mr, int n, int b) {
    int i,j,k,iam,nprocs;
    double *aux;

    aux = (double *) malloc(sizeof(double)*b*b);
    #pragma omp parallel private(iam)
    {
        nprocs=omp_get_num_threads();
        iam=omp_get_thread_num();
        
        #pragma omp for private(i,j) schedule(dynamic,1)
        for (i = 0; i < n; i += b)
        {       
                // #pragma omp for private(j) schedule(dynamic,1)
                for (j = 0; j < n; j +=b)
                {
                    reset(aux,b);   // pone a 0 la matriz auxiliar
                    for (k = 0; k < n/b; k++) // suma las tres matrices producto de la multiplicacion de bloques correspondientes
                    {
                        mm(&m1[i*n+k*b], &m2[k*b*n+j], aux, b, n);  //multiplica dos bloques

                        // //////////////DEBUG/////////////////
                        // printf("m1sub\n");
                        // escribir_submatriz(&m1[i*n+k*b],n,b);
                        // printf("m2sub\n");
                        // escribir_submatriz(&m2[k*b*n+j],n,b);
                        // printf("aux_acum\n");
                        // escribir(aux, b);
                        // ////////////////////////////////////
                    }
                    insertar_submatriz(&mr[i*n+j], aux, n, b);
                }
            
        }
    }
    free(aux);
}

int main(int argc,char *argv[]) {

    double *m1, *m2, *mr;       // matrices 1 y 2 a multiplicar, y mr es la matriz donde se guardará el resultado
    int n = DEFAULT_TAM;        // tam de la matriz cuadrada
    int b = DEFAULT_BLO;        // tam del bloque 
    int nthreads;               // # de hilos
    double start,fin;           // para el calculo de tiempos

    if (argc != 4 && argc != 1)
    {
        printf("\n\nUSO %s <dim_mat_n> <tam_blo_b> <num_hil_t>\n\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    if (argc == 4)              // si se han pasado valores en linea de comandos, se comprueban
    {
        n = atoi(argv[1]);
        b = atoi(argv[2]);
        nthreads = atoi(argv[3]);
        omp_set_num_threads(nthreads);  // se establece en número de threads
        if (n%b!=0)
        {
            printf("[ERROR] Los datos introducidos son erroneos.\nn debe ser multiplo de b.\n\n");
            exit(EXIT_FAILURE); 
        }
    }

    // omp_set_nested(1);  // para poder anidar paralelizaciones    <- NO USAR, EMPEORA MUCHO LOS RESULTADOS (SUPONGO QUE PORQUE GENERA DEMASIADOS HILOS DE EJECUCION)
    
    // inicializar matrices
    m1 = (double *) malloc(sizeof(double)*n*n); //matriz 1
    m2 = (double *) malloc(sizeof(double)*n*n); //matriz 2
    mr = (double *) malloc(sizeof(double)*n*n); //matriz resultado
    
    initializealea(m1, n*n);
    initializealea(m2, n*n);
    reset(mr, n); 

    // escribir(m1, n); // debug
    // escribir(m2, n); // debug

    start=omp_get_wtime();
    mulmatcuablo(m1, m2, mr, n, b);      //MULTIPLICACION
    // escribir(mr, n); // debug
    fin=omp_get_wtime();
    imprimir_resultados(fin-start, b, n, nthreads);

    // liberar memoria
    free(m1);
    free(m2);
    free(mr);
}
