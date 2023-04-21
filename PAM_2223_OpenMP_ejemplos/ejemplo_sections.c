#include <omp.h>
#include <stdio.h>
#include <unistd.h>

int main()  
{                                 
	int iam, np, i;
	     
	#pragma omp parallel sections private(iam, np,i)
  	{
		//printf("Hello from thread %d out of %d \n",iam,np);
        #pragma omp section
        {
			#if defined (_OPENMP) 
				np = omp_get_num_threads(); 
				iam = omp_get_thread_num();
			#endif
            printf("Soy el thread %d, actuando en solitario dentro de la seccion primera \n",iam);
            sleep(1);
        }//section
		#pragma omp section
		{
			#if defined (_OPENMP) 
				np = omp_get_num_threads(); 
				iam = omp_get_thread_num();
			#endif
			printf("Soy el thread %d, actuando en solitario dentro de la seccion segunda \n",iam);
			sleep(1);
		}//section
		#pragma omp section
		{
			#if defined (_OPENMP) 
				np = omp_get_num_threads(); 
				iam = omp_get_thread_num();
			#endif
			printf("Soy el thread %d, actuando en solitario dentro de la seccion tercera \n",iam);
			sleep(1);
		}//section
		#pragma omp section
		{
			#if defined (_OPENMP) 
				np = omp_get_num_threads(); 
				iam = omp_get_thread_num();
			#endif
			printf("Soy el thread %d, actuando en solitario dentro de la seccion cuarta \n",iam);
			sleep(1);
		}//section
    	      	
    }//parallel sections
}//main

