#include <omp.h>


#include <stdio.h>
int main()  
{                                 
        int x=1234;
	int iam, np, i,j;     
	#pragma omp parallel private(iam, np,i,x)
    	{
		#if defined (_OPENMP) 
      		  np = omp_get_num_threads(); 
      		  iam = omp_get_thread_num();
                #endif
		//printf("Hello from thread %d out of %d \n",iam,np);

         
                  printf("Soy el thread %d, antes de actualizar, con x=%d \n",iam,x);  
                  x=iam*1111;
                  printf("\t\tSoy el thread %d, despues de actualizar, con x=%d \n",iam,x);
                
                
        }//parallel
        
        printf("\n Despues de pragma parallel x=%d \n\n",x);
}//main

