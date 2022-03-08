#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// ejemplo_fibo2.c

#define COTA 20

///////////////////////////////////////////////////
long comp_fib_numbers(int n)
{
    long fnm1, fnm2, fn;
    int	iam;
  
    iam=omp_get_thread_num();
    printf(" Soy el thread %d, resolviendo fibo(%d)\n",iam,n);
  
    if ( n == 0 || n == 1 ) 
        return(1);
    /////////////////////////
    if (n < COTA)
    {
        fnm1 = comp_fib_numbers(n-1);
        fnm2 = comp_fib_numbers(n-2);
        fn = fnm1 + fnm2;
        return fn;
    }
    // Cuando una funcion de calculo de fibonacci sea menor que la cota dada,
    // no pasará de este punto de la función, por lo tanto no creará más tareas
    // y será el hilo que había cogido esta tarea el encargado de terminar el
    // calculo y las llamadas recursivas de esta funcion fibonacci en concreto
    /////////////////////////

    #pragma omp task shared(fnm1)
    {
            fnm1 = comp_fib_numbers(n-1);
    }
                            
    #pragma omp task shared(fnm2)
    {
        fnm2 = comp_fib_numbers(n-2);
    }
                            
    #pragma omp taskwait

  
    fn = fnm1 + fnm2;
                                
    return(fn);
}
                                
///////////////////////////////////////////////////
int main(int argc,char *argv[])
{
    int n=8;
    int nthreads;
    long result;
    int	iam;
  

    if (argc>1) n=atoi(argv[1]);
  
    omp_set_num_threads(4);
  
    #pragma omp parallel shared(nthreads)
    {  
        iam=omp_get_thread_num();
        printf("Soy el thread %d\n",iam);
    
        #pragma omp single nowait
        {
            result = comp_fib_numbers(n);
        } // End of single
    } // End of parallel region
          
          
    printf("Serie de Fibonacci para n=%d. El resultado es %ld \n",n,result);        
    return 0;

}