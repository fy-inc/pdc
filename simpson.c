
simpson.c

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define approx_val 2.19328059
#define N 32 /* Number of intervals in each processor */

double integrate_f(double); /* Integral function */
double simpson(int, double, double, double);

int main(int argc, char *argv[])
{
    int Procs;   /* Number of processors */
    int my_rank; /* Processor number */
    double total;
    double exact_val_of_Pi, pi, y, processor_output_share[8], x1, x2, l, sum;
    int i;
    MPI_Status status;

    /* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used. */
    MPI_Comm_size(MPI_COMM_WORLD, &Procs);

    /* Each processor computes its interval */
    x1 = ((double)my_rank) / ((double)Procs);
    x2 = ((double)(my_rank + 1)) / ((double)Procs);

    /* l is the same for all processes. */
    l = 1.0 / ((double)(2 * N * Procs));
    sum = 0.0;
    for (i = 1; i < N; i++)
    {
        y = x1 + (x2 - x1) * ((double)i) / ((double)N);

        /* call Simpson's rule  */
        sum = (double)simpson(i, y, l, sum);
    }

    /* Include the endpoints of the intervals */
    sum += (integrate_f(x1) + integrate_f(x2)) / 2.0;
    total = sum;

    /* Add up the integrals calculated by each process. */
    if (my_rank == 0)
    {
        processor_output_share[0] = total;

        /* source = i, tag = 0 */
        for (i = 1; i < Procs; i++)
            MPI_Recv(&(processor_output_share[i]), 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
    }
    else
    {
        /* dest = 0, tag = 0 */
        MPI_Send(&total, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    /* Add up the value of Pi and print the result.  */
    if (my_rank == 0)
    {
        for (i = 0; i < Procs; i++)
            pi += processor_output_share[i];
        pi *= 2.0 * l / 3.0;
        printf("-------------------------------------------------\n");
        printf("The computed Pi of the integral for %d grid points is  %25.16e \n",
               (N * Procs), pi);

        /* This is directly derived from the integeration of the formula. See 
  the report. */
#if 1
        exact_val_of_Pi = 4.0 * atan(1.0);
#endif

#if 0
      exact_val_of_Pi = 4.0 * log(approx_val);
#endif
        printf("The error or the discrepancy between exact and computed value of Pi : %25.16e\n",
               fabs(pi - exact_val_of_Pi));
        printf("-------------------------------------------------\n");
    }

    MPI_Finalize();
}

double integrate_f(double x)
{
    /* compute and return value */
    return 4.0 / (1.0 + x * x);
}

double simpson(int i, double y, double l, double sum)
{
    /* store result in sum */
    sum += integrate_f(y);
    sum += 2.0 * integrate_f(y - l);
    if (i == (N - 1))
        sum += 2.0 * integrate_f(y + l);
    return sum;
} /* simpson */

/*
#include "mpi.h"
#include<stdlib.h>
#include<unistd.h>
#include<stdio.h>
#include<math.h>

int n;
int N;
int epp;
float f(float x,int i)
{
    if(i==0)
        return 1.0/(x*x+1);
    if(i==n)
        return 1.0/(x*x+1);
    if(i%2)
        return 4.0/(x*x+1);
    return 2.0/(x*x+1);
}

int main(int argc,char *argv[])
{
    int pid,n_er,np;
    float a,b;
    //scanf("%d%f%f",&n,&a,&b);
    n = 6;
    a = 0;
    b = 1;
    N = n+1;
    float h = (b-a)/(float)n;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    epp = N/np;
    if(pid==0)
    {
        //printf("%d\n",pid);
        //printf("%d %d %d\n",epp,N, np);
        int index,i;
        if(np>1)
        {
            for(i=1;i<np-1;i++)
                MPI_Send(&epp,1,MPI_INT,i,0,MPI_COMM_WORLD);
            //printf("%d",i);   
            index = i*epp;
            int e_left = N - index;
            MPI_Send(&e_left,1,MPI_INT,i,0,MPI_COMM_WORLD);
        }
        float sum = 0;
        for(i=0;i<epp;i++)
        {
            float val = f(a+((float)i)*h,i);
            sum = sum + val;
            //printf("%f %d",val,i);
        }
        for(i=1;i<np;i++)
        {
            float temp;
            MPI_Recv(&temp,1,MPI_FLOAT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
            sum = sum + temp;
        }
        float ans = (h*sum)/3.0;
        printf("the value is: %f",ans);
    }
    else
    {
        //printf("%d\n",pid);
        MPI_Recv(&n_er,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        float partial_sum = 0;
        int i;
        for(i=0;i<n_er;i++)
        {
            float index = pid*epp + i;
            //printf("%f %d %d\n",index,epp,i);
            float val = f(a+index*(h),index);
            partial_sum +=val;
            //printf("%f %f \n",val,index);
            //printf("%f",partial_sum);
        }
        MPI_Send(&partial_sum,1,MPI_FLOAT,0,0,MPI_COMM_WORLD);    
    }
    MPI_Finalize();
    return 0;
}
*/



