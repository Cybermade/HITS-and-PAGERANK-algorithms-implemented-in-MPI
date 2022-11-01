/*
 *
 * mpicc -o hits_algorithm_mpi hits_algorithm_mpi.c
 * mpirun -np 8 ./hits_algorithm_mpi
 *
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void read_ints(const char *file_name, double *graph, int nb_nodes)
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

    int i = 0;
    const int nb_nodes = 8;
    double *buff_hub;
    double *buff_autority;
    double *graph;
    double *hub;
    double *autority;
    int nb_iterations = 10;

    int *sendcounts;               // array describing how many elements to send to each process
    int *displs;                   // array describing the displacements where each segment begins
    int rem = nb_nodes % numtasks; // elements remaining after division among processes
    int sum = 0;

    graph = calloc(nb_nodes * nb_nodes, sizeof(double));
    hub = calloc(nb_nodes, sizeof(double));
    autority = calloc(nb_nodes, sizeof(double));
    sendcounts = calloc(numtasks, sizeof(int));
    displs = calloc(numtasks, sizeof(int));

    for (int i = 0; i < nb_nodes; i++)
    {

        autority[i] = 1;
        hub[i] = 1;
    }

    if (rank == 0)
    {
        read_ints("test2.txt", graph, nb_nodes);
        for (int i = 0; i < nb_nodes; i++)
        {
            hub[i] = 1;
            autority[i] = 1;
        }
    }

    MPI_Bcast(graph, nb_nodes * nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("i am %d\n", rank);
    for (int i = 0; i < nb_nodes; i++)
    {
        for (int j = 0; j < nb_nodes; j++)
        {
            printf("%.5f ", graph[i * nb_nodes + j]);
        }
        printf("\n");
    }
    printf("\n\n\n");

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

    buff_hub = calloc(sendcounts[rank], sizeof(int));
    buff_autority = calloc(sendcounts[rank], sizeof(int));
    
    MPI_Bcast(hub, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(autority, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double norm = 0;

    int begin = displs[rank];
    for (int k = 0; k < nb_iterations; k++)
    {
        for (int i = 0; i < sendcounts[rank]; i++)
        {
            buff_autority[i] = 0;
            for (int j = 0; j < nb_nodes; j++)
            {
                if (graph[j * nb_nodes + i + begin])
                {
                    buff_autority[i] += hub[j];
                }
            }
        }
        MPI_Gatherv(buff_autority, sendcounts[rank], MPI_DOUBLE, autority, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            norm = 0;
            for (int i = 0; i < nb_nodes; i++)
            {
                norm += pow(autority[i], 2);
            }
            norm = sqrt(norm);
               
        }
        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for(int i = 0; i<sendcounts[rank]; i++){
            buff_autority[i] = buff_autority[i]/norm;
        }
        MPI_Gatherv(buff_autority, sendcounts[rank], MPI_DOUBLE, autority, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(autority, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = 0; i < sendcounts[rank]; i++)
        {
            buff_hub[i] = 0;
            for (int j = 0; j < nb_nodes; j++)
            {
                if (graph[(i + begin) * nb_nodes + j])
                {
                    buff_hub[i] += autority[j];
                }
            }
        }
        
        MPI_Gatherv(buff_hub, sendcounts[rank], MPI_DOUBLE, hub, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            norm = 0;
            for (int i = 0; i < nb_nodes; i++)
            {
                norm += pow(hub[i], 2);
            }
            norm = sqrt(norm);
            
        }
        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for(int i = 0; i<sendcounts[rank]; i++){
            buff_hub[i] = buff_hub[i]/norm;
        }
        MPI_Gatherv(buff_hub, sendcounts[rank], MPI_DOUBLE, hub, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(hub, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (rank == 0)
        {
            printf("Iteration NB : %d\nHub :\t\tAutority :\n",k);
            for (int i = 0; i < nb_nodes; i++)
            {
                printf("%d -> %.5f , %d -> %.5f\n", i, hub[i], i, autority[i]);
            }
            printf("\n");
        }
    }
    

    free(graph);
    free(hub);
    free(autority);
    free(sendcounts);
    free(displs);
    free(buff_hub);
    free(buff_autority);

    // done with MPI
    MPI_Finalize();
}