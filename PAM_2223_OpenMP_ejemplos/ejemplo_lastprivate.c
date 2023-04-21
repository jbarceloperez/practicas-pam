#include <omp.h>


#include <stdio.h>
int main()  
{                                 
        int x=9999;
	int iam =0, np = 1, i=0,j=0;     
	
	
	printf("\n Antes de pragma parallel x=%d \n\n",x);
	
	#pragma omp parallel private(iam, np,i) 
    	{
		#if defined (_OPENMP) 
      		  np = omp_get_num_threads(); 
      		  iam = omp_get_thread_num();
                #endif
		//printf("Hello from thread %d out of %d \n",iam,np);

         
                  printf("\t Soy el thread %d, antes de entrar en bucle, con x=%d \n",iam,x);  // mensajito para que le de tiempo al master a poner x=999             
                  
                  #pragma omp for lastprivate(x) schedule(dynamic)
                  for(i=0;i<11;i++)
                  {
                      printf("\t\t\t Soy el thread %d, antes de actualizar en for, i=%d x=%d\n",iam,i,x);
                      x=iam*i;
                      printf("\t\t\t\t Soy el thread %d, actualizando en for, i=%d x=%d\n",iam,i,x);
                  }  
                  printf("\t Soy el thread %d, despues de salir de bucle, con x=%d \n",iam,x); 
                  //todos se quedan con el valor que ponga el thread que haga la iteracion i=10
                
                
        }//parallel
        
        printf("\n Despues de pragma parallel x=%d \n\n",x);
}//main

