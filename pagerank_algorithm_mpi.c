/*
 *
 * mpicc -o hits_algorithm_mpi hits_algorithm_mpi.c -lm
 * mpirun -np 8 ./hits_algorithm_mpi
 *
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void read_ints(const char *file_name, int *graph, int nb_nodes)
{
    FILE *file = fopen(file_name, "r");
    int i = 0;
    int j = 0;

    while (!feof(file))
    {

        fscanf(file, "%d", &i);
        fscanf(file, "%d", &j);
        graph[(int)i * nb_nodes + j] = 1;
    }
    fclose(file);
}
void update_pagerank(double *pagerank, double *pagerank_old, int nb_nodes)
{
    for(int i = 0; i < nb_nodes; i++)
    {
        pagerank_old[i] = pagerank[i];
        pagerank[i] = 0;
    }
}
int number_of_neighbors(int *graph, int node, int nb_nodes)
{
    int sum = 0;
    for (int i = 0; i < nb_nodes; i++)
    {   
        if(graph[node * nb_nodes + i])
        {
            sum++;
        }
        
        
    }
    return sum;
}
int main(int argc, char *argv[])
{
    int numtasks, rank, len, rc;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    char message[] = "Hello!";
    int number[2];
    // initialize MPI
    MPI_Init(&argc, &argv);

    // get number of tasks
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    // get my rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // this one is obvious
    MPI_Get_processor_name(hostname, &len);
    // printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks ,rank ,hostname);

    double *pagerank = NULL;
    int *graph = NULL;
    int nb_nodes = 4;
    
    int *sendcounts;               // array describing how many elements to send to each process
    int *displs;                   // array describing the displacements where each segment begins
    int rem = nb_nodes % numtasks; // elements remaining after division among processes
    int sum = 0;
    
    double *pagerank_old = NULL;
    int nb_iterations = 5;
    graph = calloc(nb_nodes *nb_nodes , sizeof (int));
    
    pagerank_old = calloc(nb_nodes, sizeof (double));
    sendcounts = calloc(numtasks, sizeof (int));
    displs = calloc(numtasks, sizeof (int));
    double *buff_pagerank = NULL;

    


    

    if (rank == 0)
    {
        read_ints("test3.txt", graph, nb_nodes);
        
        pagerank = calloc(nb_nodes, sizeof (double));
        for(int i=0;i<nb_nodes;i++)
        {
            pagerank[i] = 1.0/nb_nodes;
        }
        

    }
    
    MPI_Bcast(graph, nb_nodes * nb_nodes, MPI_INT, 0, MPI_COMM_WORLD);

    // printf("i am %d\n", rank);
    // for (int i = 0; i < nb_nodes; i++)
    // {
    //     for (int j = 0; j < nb_nodes; j++)
    //     {
    //         printf("%d ", graph[i * nb_nodes + j]);
    //     }
    //     printf("\n");
    // }


    for (int i = 0; i < numtasks; i++)
    {
        sendcounts[i] = nb_nodes / numtasks;
        if (rem > 0)
        {
            sendcounts[i]++;
            rem--;
        }
        displs[i] = sum;
        sum += sendcounts[i];
    }
    buff_pagerank = calloc(sendcounts[rank], sizeof (double));

    int begin = displs[rank];
    
    for (int k =0;k<=nb_iterations;k++)
    {  if(rank==0){
        printf("iteration %d\n",k);
        for(int i=0;i<nb_nodes;i++)
        {
            printf("%d -> %.5f\n",i,pagerank[i]);
        }
        double sum = 0;
        for (int i = 0; i < nb_nodes; i++)
        {
            sum+= pagerank[i];
            
        }
        printf("Sum : %.5f (Should always be 1)\n\n", sum);
        printf("\n\n\n");


        update_pagerank(pagerank, pagerank_old, nb_nodes);
    }
        MPI_Bcast(pagerank_old, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = 0; i < sendcounts[rank]; i++)
        {
            buff_pagerank[i] = 0;
            //printf("%d : \n",i);
            for (int j = 0; j < nb_nodes; j++)
            {   
                if (graph[j * nb_nodes + i + begin])
                {   
                    //printf("%d %d\n",i+begin,j);
                    buff_pagerank[i] += (pagerank_old[j]/(double)number_of_neighbors(graph, j, nb_nodes));
                }
            }
            
        }
        MPI_Gatherv(buff_pagerank, sendcounts[rank], MPI_DOUBLE, pagerank, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        
        
    }
    
    


    free(graph);
    if(rank==0)free(pagerank);
    free(pagerank_old);
    // done with MPI
    MPI_Finalize();
    
    
    

    
   


}