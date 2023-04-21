#include <omp.h>
#include <stdio.h>

int main()
{
        int data;
        int flag=0;
        int iam;

        #pragma omp parallel sections num_threads(2)  private(iam)
        {
                #pragma omp section
                {
                        iam=omp_get_thread_num();
                        printf("Thread %d, esperando dato: \n",iam);
                        scanf("%d",&data);

                        #pragma omp flush(data)                        
			       flag = 1;
                        #pragma omp flush(flag)
			}


		#pragma omp section
                {
                        iam=omp_get_thread_num();
                        printf("ANTES: Thread %d, dato a procesar = %d \n",iam,data);

                        while (!flag)
                        {
                                #pragma omp flush(flag)
                        }
                        #pragma omp flush(data)

                        printf("DESPUES: Thread %d: dato a procesar = %d \n", iam,data);

                        // do more work with data
                }
        } //parallel
}


