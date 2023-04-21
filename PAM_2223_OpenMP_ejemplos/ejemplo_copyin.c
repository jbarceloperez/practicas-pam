#include <omp.h>


#include <stdio.h>


int x;
#pragma omp threadprivate(x) // privada a cada thread --> no se copia el valor del master 

int main()  
{                                 
	int iam =0, np = 1, i=0,j=0;     
	
	x=9999; // lo ponemos en el master
	
	#pragma omp parallel private(iam, np,i) copyin(x) //sin el copyin, en cada thread salvo el master, el valor es x=0
    	{
		#if defined (_OPENMP) 
      		  np = omp_get_num_threads(); 
      		  iam = omp_get_thread_num();
                #endif
		//printf("Hello from thread %d out of %d \n",iam,np);

		  sleep(1);
                  printf("Soy el thread %d, antes de actualizar, con x=%d \n",iam,x);  // mensajito para que le de tiempo al master a poner x=999             
                  x=iam;
                  printf("\t\tSoy el thread %d, despues de que el master actualice, con x=%d \n",iam,x); 
                  
                
        }//parallel
        
        printf("\n Despues de pragma parallel x=%d \n\n",x);
}//main

