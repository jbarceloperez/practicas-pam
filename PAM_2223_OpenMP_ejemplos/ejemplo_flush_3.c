#include <omp.h>
#include <stdio.h>

int main()
{
    int data1, data2, suma;
    int flag1=0, flag2=0;
    int iam;

    #pragma omp parallel sections num_threads(3)  private(iam)
    {
        #pragma omp section
        {
            iam=omp_get_thread_num();
            printf("Thread %d, esperando primer dato: \n",iam);
            scanf("%d",&data1);
            printf("Thread %d, esperando segundo dato: \n",iam);
            scanf("%d",&data2);

            #pragma omp flush(data1, data2)                        
            flag1 = 1;
            #pragma omp flush(flag1)
        } //section

        #pragma omp section
        {
            iam=omp_get_thread_num();
            printf("ANTES: Thread %d, esperando datos...\n",iam);

            while (!flag1)
            {
                #pragma omp flush(flag1)
            }
            #pragma omp flush(data1, data2)

            printf("DESPUES: Thread %d: sumando %d + %d \n", iam,data1,data2);

            suma = data1 + data2;
            #pragma omp flush(suma)
            flag2 = 1;
            #pragma omp flush(flag2)
        } //section

        #pragma omp section
        {
            iam=omp_get_thread_num();
            printf("ANTES: Thread %d, esperando suma...\n",iam);

            while (!flag2)
            {
                #pragma omp flush(flag2)
            }
            #pragma omp flush(suma)

            printf("DESPUES: Thread %d: la suma es %d\n", iam, suma);

        } //section

    } //parallel
}