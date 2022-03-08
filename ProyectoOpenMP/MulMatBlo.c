// SecMulMatBlo.c
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h> 

#define DEFAULT_TAM 16
#define DEFAULT_BLO 4

// imprime la matriz, para debuggear
void escribir(double *mr, int m, int n) {
  int i,j;
  for (i=0;i<m;i++) {
    for(j=0;j<n;j++)
      printf("%.4lf\t",mr[i*n+j]);
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
void imprimir_resultados(double tiempo, int b, int m, int n, int nthreads) {
    double mflops; 

    if (tiempo==0.)
    {
      printf("No hay suficiente precision\n");
    }
    else 
    {
      mflops=((2.*m*n)/tiempo)/1000000.;
      printf("\n  Numero de hilos: %d\n  Tam de matriz: %d x %d\n  Tam de bloque: %d\n    segundos: %.6lf, Mflops: %.6lf\n",nthreads,m,n,b,tiempo,mflops);
    }
}

/**
 *  Pone la matriz entera a 0
 *  @param b Tam de la matriz
 */
void reset(double *m, int b) {
    int i;
    for (i = 0; i < b; i++)
        m[i] = 0.;
}

/** 
 * Inicializa un BLOQUE con valores aleatorios
 * @param t Tam de la matriz 
 */ 
void initializealea(double *m, int t) {     
    int i;
    for (i = 0; i < t; i++)
        m[i] = (double)rand()/RAND_MAX;
}

// multiplica dos submatrices 
void mm(double *m1, double *m2, double *mr, int t, int m, int n, int k) {
    int i,j,h;
    double suma;
 
    #pragma omp parallel for private(i,j,suma) schedule(dynamic, 1)
    for (i = 0; i < t; i++) {
        for(j=0;j<t;j++){
            suma=mr[i*t+j]; // se obtiene el valor guardado en aux para seguir sumándoselo
            for(h=0;h<t;h++)
                suma += m1[i*k+h] * m2[j+h*n];  
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

void mulmatcuablosec(double *m1, double *m2, double *mr, int m, int n, int k, int b) {
    int i,j,h,iam,nprocs;
    double *aux;

    aux = (double *) malloc(sizeof(double)*b*b);
    #pragma omp parallel private(iam)
    {
        nprocs=omp_get_num_threads();
        iam=omp_get_thread_num();

        #pragma omp for private(i,j) schedule(dynamic,1)
        for (i = 0; i < m; i += b)
        {
            for (j = 0; j < n; j +=b)
            {
                reset(aux,b*b);   // pone a 0 la matriz auxiliar
                for (h = 0; h < k/b; h++) // suma las matrices producto de la multiplicacion de bloques correspondientes
                {
                    mm(&m1[i*k+h*b], &m2[h*b*n+j], aux, b, m, n, k);  //multiplica dos bloques

                    // //////////////DEBUG/////////////////
                    // printf("m1sub\n");
                    // escribir_submatriz(&m1[i*n+h*b],n,b);
                    // printf("m2sub\n");
                    // escribir_submatriz(&m2[h*b*n+j],n,b);
                    // printf("aux_acum\n");
                    // escribir(aux, b, b);
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
    int m = DEFAULT_TAM;        // tam m de la matriz resultado (# de filas)
    int n = DEFAULT_TAM;        // tam n de la matriz resultado (# de columnas)
    int k = DEFAULT_TAM;        // tam k de las matrices de entrada (# de columnas de la primera y # de filas de la segunda) 
    int b = DEFAULT_BLO;        // tam del bloque
    int nthreads;               // # de hilos
    double start,fin;           // para el calculo de tiempos

    if (argc != 6 && argc != 1)
    {
        printf("\n\nUSO %s <dim_mat_m> <dim_mat_n> <dim_mat_k> <tam_blo_b> <num_hil_t>\n\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    if (argc == 6)              // si se han pasado valores en linea de comandos, se comprueban
    {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        k = atoi(argv[3]);
        b = atoi(argv[4]);
        nthreads = atoi(argv[5]);
        omp_set_num_threads(nthreads);  // se establece en número de threads
        if (n%b!=0 || m%b!=0 || k%b!=0)
        {
            printf("[ERROR] Los datos introducidos son erroneos.\nLas dimensiones de las matrices deben ser todas multiplos del tam de bloque.\n\n");
            exit(EXIT_FAILURE); 
        }
    }

    // inicializar matrices
    m1 = (double *) malloc(sizeof(double)*m*k); //matriz 1
    m2 = (double *) malloc(sizeof(double)*k*n); //matriz 2
    mr = (double *) malloc(sizeof(double)*m*n); //matriz resultado
    
    initializealea(m1, m*k);
    initializealea(m2, k*n);
    reset(mr, m*n); 

    // escribir(m1, m, k); // debug
    // escribir(m2, k, n); // debug

    start=omp_get_wtime();
    mulmatcuablosec(m1, m2, mr, m, n, k, b);      //MULTIPLICACION
    // escribir(mr, m, n); // debug
    fin=omp_get_wtime();
    imprimir_resultados(fin-start, b, m, n, nthreads);

    // liberar memoria
    free(m1);
    free(m2);
    free(mr);
}
