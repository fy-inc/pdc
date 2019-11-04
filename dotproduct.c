
dotproduct.c

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define VECLEN 100

int main(int argc, char *argv[])
{
    int i, myid, numprocs, len = VECLEN;
    double *x, *y;
    double mysum, allsum;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    if (myid == 0)
        printf("Starting omp_dotprod_mpi. Using %d tasks...\n", numprocs);


    x = (double *)malloc(len * sizeof(double));
    y = (double *)malloc(len * sizeof(double));

    for (i = 0; i < len; i++)
    {
        x[i] = 1.0;
        y[i] = x[i];
    }

    mysum = 0.0;
    for (i = 0; i < len; i++)
    {
        mysum += x[i] * y[i];
    }

    printf("Task %d partial sum = %f\n", myid, mysum);

    MPI_Reduce(&mysum, &allsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid == 0)
        printf("Done. MPI version: global sum  =  %f \n", allsum);

    free(x);
    free(y);
    MPI_Finalize();
}
/*
#include "mpi.h"
#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>

#define n 10

int arr[] = {1,2,3,4,5,6,7,8,9,10};
int temp_arr[1000];
int brr[] = {1,2,3,4,5,6,7,8,9,10};
int temp_brr[1000];

int main(int argc,char *argv[])
{
    int pid,np,n_er;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    printf("%d\n",pid);
    if(pid==0)
    {
        int index,i;
        int epp = n/np;
        if(np>1)
        {
            for(i=1;i<np-1;i++)
            {
                index = i*epp;
                MPI_Send(&epp,1,MPI_INT,i,0,MPI_COMM_WORLD);
                MPI_Send(&arr[index],epp,MPI_INT,i,0,MPI_COMM_WORLD);
                MPI_Send(&brr[index],epp,MPI_INT,i,0,MPI_COMM_WORLD);
            }
            index = i*epp;
            int elements_left = n - index;
            MPI_Send(&elements_left,1,MPI_INT,i,0,MPI_COMM_WORLD);
            MPI_Send(&arr[index],elements_left,MPI_INT,i,0,MPI_COMM_WORLD);
            MPI_Send(&brr[index],elements_left,MPI_INT,i,0,MPI_COMM_WORLD);
        }
        int sum = 0;
        for(i=0;i<epp;i++)
            sum += arr[i]*brr[i];
        for(i=1;i<np;i++)
        {
            int temp;
            MPI_Recv(&temp,1,MPI_INT,MPI_ANY_SOURCE,0, MPI_COMM_WORLD,&status);
            sum = sum + temp;
        }
        printf("sum is: %d",sum);
    }
    else
    {
        MPI_Recv(&n_er,1,MPI_INT, 0,0,MPI_COMM_WORLD,&status);
        MPI_Recv(&temp_arr,n_er,MPI_INT, 0,0,MPI_COMM_WORLD,&status);
        MPI_Recv(&temp_brr,n_er,MPI_INT, 0,0,MPI_COMM_WORLD,&status);
        int partial_sum = 0;
        for(int i=0;i<n_er;i++)
            partial_sum+= temp_arr[i]*temp_brr[i];
        MPI_Send(&partial_sum,1,MPI_INT,0,0,MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
*/



