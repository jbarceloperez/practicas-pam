#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

//ejemplo_copyin.c

int x;
#pragma omp threadprivate(x) // privada a cada thread --> no se copia el valor del master 


int main()  
{                                 
	int iam,np,i,j;     
	
	x=9999; // lo ponemos en el master
	
	#pragma omp parallel private(iam, np,i) copyin(x) 
	//si no se pusiera el copyin(x), en cada thread salvo el master, el valor sería x=0 por ser una variable definida como prevada de cda thread de manera general
    	{
		#if defined (_OPENMP) 
      		  np = omp_get_num_threads(); 
      		  iam = omp_get_thread_num();
                #endif
		//printf("Hello from thread %d out of %d \n",iam,np);

                 printf("Soy el thread %d, antes de actualizar, con x=%d \n",iam,x); 
                 x=iam+1000;
                 printf("\t\tSoy el thread %d, despues de actualizar, con x=%d \n",iam,x); 
                  
                
        }//parallel
        
        printf("\n Despues de pragma parallel x=%d \n\n",x);
}//main

