/*
 *
 * @author : Cybermade
 * mpicc -o hits_algorithm_mpi hits_algorithm_mpi.c -lm
 * mpirun -np <nb_of_tasks> ./hits_algorithm_mpi <input_file> <nb_nodes> <nb_iterations>
 *
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* functions prototypes */
void read_ints(const char *, int *, int);      /* read graph matrix from from file */

int main(int argc, char *argv[])
{
    /* MPI variables */
    int numtasks, rank, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    
    /* MPI initialization */
    MPI_Init(&argc, &argv);

    /* Get the number of tasks */
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    /* Get the rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* this one is obvious */
    MPI_Get_processor_name(hostname, &len);
   
    /* Variables */
    
    char *nbnodes = argv[2];
    int nb_nodes = atoi(nbnodes); /* number of nodes */
    

    int *graph = NULL; /* graph */
    double *hub = NULL;     /* hub */
    double *autority = NULL;     /* autority */
    double *buff_hub = NULL;           /* buffer for hub vector (used for MPI communication) */
    double *buff_autority = NULL;      /* buffer for autority vector (used for MPI communication) */
    
    char *nbiterations = argv[3];
    int nb_iterations = atoi(nbiterations); /* number of iterations to run the algorithm */

    int *sendcounts;               /* array describing how many elements to send to each process */
    int *displs;                   /* array describing the displacements where each segment begins */
    int rem = nb_nodes % numtasks; /* elements remaining after division among processes */
    int sum = 0;
    int begin = 0;                 /* begin index for each process */
    double norm = 0;               /* normalization variable */

    /* memory allocation */
    graph = calloc(nb_nodes * nb_nodes, sizeof(int));
    hub = calloc(nb_nodes, sizeof(double));
    autority = calloc(nb_nodes, sizeof(double));
    sendcounts = calloc(numtasks, sizeof(int));
    displs = calloc(numtasks, sizeof(int));

    /* read graph from file (we assume that the root process has the graph and it is broadcasted to all processes) */
    if (rank == 0)
    {
        /* read graph from file */
        read_ints(argv[1], graph, nb_nodes); /* read graph from file */
        /* initialize hub and autority vectors */
        for (int i = 0; i < nb_nodes; i++)
        {
            hub[i] = 1;
            autority[i] = 1;
        }
        /* calculate sendcounts and displs */
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
        /* create or clear the output file */
        fclose(fopen("result_hits_mpi.txt", "w"));
    }
    /* broadcast graph, hub, autority to all processes */
    MPI_Bcast(graph, nb_nodes * nb_nodes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(hub, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(autority, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    
    

    /* broadcast sendcounts and displs to all processes */
    MPI_Bcast(sendcounts, numtasks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, numtasks, MPI_INT, 0, MPI_COMM_WORLD);

    /* allocate memory for buffers */
    buff_hub = calloc(sendcounts[rank], sizeof(double));
    buff_autority = calloc(sendcounts[rank], sizeof(double));

    /* calculate begin index for each process */
    begin = displs[rank];

    
    
    /* begin computation */
    for (int k = 0; k < nb_iterations; k++)
    {   
        
        
        /* calculate autority vector */
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
        /* gather autority vector */
        MPI_Gatherv(buff_autority, sendcounts[rank], MPI_DOUBLE, autority, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        /* root calculate sum of autority vector */
        if (rank == 0)
        {
            norm = 0;
            for (int i = 0; i < nb_nodes; i++)
            {
                norm += pow(autority[i], 2);
            }
            norm = sqrt(norm);
               
        }
        /* broadcast normalize variable to all processes */
        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /* each process normalize its part of the autority vector */
        for(int i = 0; i<sendcounts[rank]; i++){
            buff_autority[i] = buff_autority[i]/norm;
        }
        /* print hub and autority vectors */
        
        /* gather autority vector */
        MPI_Gatherv(buff_autority, sendcounts[rank], MPI_DOUBLE, autority, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /* broadcast autority vector to all processes */
        MPI_Bcast(autority, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /* calculate hub vector */
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
        /* gather hub vector */
        MPI_Gatherv(buff_hub, sendcounts[rank], MPI_DOUBLE, hub, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /* root calculate sum of hub vector */
        if (rank == 0)
        {
            norm = 0;
            for (int i = 0; i < nb_nodes; i++)
            {
                norm += pow(hub[i], 2);
            }
            norm = sqrt(norm);
            
        }
        
        /* broadcast normalize variable to all processes */
        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /* each process normalize its part of the hub vector */
        for(int i = 0; i<sendcounts[rank]; i++){
            buff_hub[i] = buff_hub[i]/norm;
        }
        /* gather hub vector */
        MPI_Gatherv(buff_hub, sendcounts[rank], MPI_DOUBLE, hub, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /* broadcast hub vector to all processes */
        MPI_Bcast(hub, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (rank == 0)
        {
            FILE *fptr;
            fptr = fopen("result_hits_mpi.txt", "a");
            printf("Iteration NB : %d\nHub :\t\tAutority :\n",k);
            fprintf(fptr, "Iteration NB : %d\nHub :\t\tAutority :\n",k);
            for (int i = 0; i < nb_nodes; i++)
            {
                printf("%d -> %.5f , %d -> %.5f\n", i, hub[i], i, autority[i]);
                fprintf(fptr, "%d -> %.5f , %d -> %.5f\n", i, hub[i], i, autority[i]);
            }
            printf("\n\n\n");
            fprintf(fptr, "\n\n\n");
            
            fclose(fptr);
        }
    }
    
    /* free memory */
    free(graph);
    graph = NULL;
    free(hub);
    hub = NULL;
    free(autority);
    autority = NULL;
    free(sendcounts);
    sendcounts = NULL;
    free(displs);
    displs = NULL;
    free(buff_hub);
    buff_hub = NULL;
    free(buff_autority);
    buff_autority = NULL;

    /* finalize MPI */
    MPI_Finalize();
}
void read_ints(const char *file_name, int *graph, int nb_nodes)
{
    if (file_name == NULL)
    {
        printf("Error: file name is NULL\n");
        exit(1);
    }
    /* open file */
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