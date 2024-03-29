/*
 *
 * @author : Cybermade
 * mpicc -Wall -o pagerank_algorithm_mpi pagerank_algorithm_mpi.c
 * mpirun -np <nb_of_tasks> ./pagerank_algorithm_mpi <input_file> <nb_nodes> <nb_iterations>
 *
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* functions prototypes */
void read_ints(const char *, int *, int);      /* read graph matrix from from file */
void update_pagerank(double *, double *, int); /* save the current pagerank before going to the next iteration */
int number_of_neighbors(int *, int, int);      /* return the number of neighbors of a node */

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
    int *graph = NULL; /* graph */
    char *nbnodes = argv[2];
    int nb_nodes = atoi(nbnodes); /* number of nodes */

    double *pagerank = NULL;     /* pagerank */
    double *pagerank_old = NULL; /* old pagerank vector (we need this since we calculate pagerank using the last iteration) */
    char *nbiterations = argv[3];
    int nb_iterations = atoi(nbiterations); /* number of iterations to run the algorithm */
    double *buff_pagerank = NULL;           /* buffer for pagerank vector (used for MPI communication) */

    int *sendcounts = NULL;        /* array describing how many elements to send to each process */
    int *displs = NULL;            /* array describing the displacements where each segment begins */
    int rem = nb_nodes % numtasks; /* elements remaining after division among processes */
    int sum = 0;                   /* sum of sendcounts */
    double dumping_factor = 0.85;  /* dumping factor for pagerank */
    int begin = 0;                 /* begin index for each process */

    /* memory allocation */
    graph = calloc(nb_nodes * nb_nodes, sizeof(int));
    pagerank_old = calloc(nb_nodes, sizeof(double));
    sendcounts = calloc(numtasks, sizeof(int));
    displs = calloc(numtasks, sizeof(int));
    

    /* read graph from file (we assume that the root process has the graph and it is broadcasted to all processes) */
    if (rank == 0)
    {
        read_ints(argv[1], graph, nb_nodes); /* read graph from file */

        /*
         *  alocate memory for pagerank, this is needed only here since only the root needs to save the current pagerank iteration
         *  (nodes use pagerank_old to calculate the new pagerank)
         */

        pagerank = calloc(nb_nodes, sizeof(double));

        /* initialize pagerank */
        for (int i = 0; i < nb_nodes; i++)
        {
            pagerank[i] = 1.0 / nb_nodes;
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
        fclose(fopen("result_pagerank_mpi.txt", "w"));
    }

    /* broadcast graph to all processes */
    MPI_Bcast(graph, nb_nodes * nb_nodes, MPI_INT, 0, MPI_COMM_WORLD);

    /* broadcast sendcounts and displs */
    MPI_Bcast(sendcounts, numtasks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, numtasks, MPI_INT, 0, MPI_COMM_WORLD);
    /* allocate memory for pagerank buffer */
    buff_pagerank = calloc(sendcounts[rank], sizeof(double));

    /* calculate begin index for each process */
    begin = displs[rank];

    /* begin computation */
    for (int k = 0; k <= nb_iterations; k++)
    {
        /* only root updates the pagerank and prints the pagerank after each iteration */
        if (rank == 0)
        {   
            FILE *fptr;
            fptr = fopen("result_pagerank_mpi.txt", "a");
            printf("iteration %d\n", k);
            fprintf(fptr, "iteration %d\n", k);
            for (int i = 0; i < nb_nodes; i++)
            {
                printf("%d -> %.5f\n", i, pagerank[i]);
                fprintf(fptr, "%d -> %.5f\n", i, pagerank[i]);
            }
            double sum = 0;
            fprintf(fptr, "\n");
            for (int i = 0; i < nb_nodes; i++)
            {
                sum += pagerank[i];
            }
            printf("Sum : %.5f (From 0 to 1, the higher the better for a web structure)\n\n", sum);
            fprintf(fptr, "Sum : %.5f (From 0 to 1, the higher the better for a web structure)\n\n", sum);
            printf("\n\n\n");
            fprintf(fptr, "\n\n\n");

            /* update pagerank */
            update_pagerank(pagerank, pagerank_old, nb_nodes);

            fclose(fptr);
        }
        /* broadcast pagerank_old to all processes */
        MPI_Bcast(pagerank_old, nb_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /* calculate pagerank for each node */
        for (int i = 0; i < sendcounts[rank]; i++)
        {
            buff_pagerank[i] = (1 - dumping_factor) / nb_nodes;

            for (int j = 0; j < nb_nodes; j++)
            {
                if (graph[j * nb_nodes + i + begin])
                {
                    buff_pagerank[i] += ((dumping_factor) * (pagerank_old[j] / (double)number_of_neighbors(graph, j, nb_nodes)));
                }
            }
        }
        /* gather pagerank from all processes */
        MPI_Gatherv(buff_pagerank, sendcounts[rank], MPI_DOUBLE, pagerank, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    if(rank ==0)
    {   
        /* simple algorithm to get the top pages */
        FILE *fptr;
        fptr = fopen("result_pagerank_mpi.txt", "a");
        double max = pagerank[0]; 
        int max_index = 0;
        int k = nb_nodes;
        printf("Top pages are -> pagerank (descending order) :\n");
        fprintf(fptr, "Top pages are -> pagerank (descending order) :\n");
        for(int i = 0;i<nb_nodes;i++)
        {   
            
            for(int j = 0; j < nb_nodes; j++)
            {
                if(pagerank[j] > max)
                {
                    max = pagerank[j];
                    max_index = j;
                }
                
            }
            
            printf("%d -> %d\n", max_index, k);
            fprintf(fptr, "%d -> %d\n", max_index, k);
            k--;
            pagerank[max_index] = -1;
            max = -1;
            
        }
        fclose(fptr);
    }
    /* free memory */
    printf("Process %d finished\n", rank);
    free(graph);
    graph = NULL;
    if (rank == 0)
        {free(pagerank);pagerank = NULL;}
    free(pagerank_old);
    pagerank_old = NULL;
    free(sendcounts);
    sendcounts = NULL;
    free(displs);
    displs = NULL;
    free(buff_pagerank);
    buff_pagerank = NULL;


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
void update_pagerank(double *pagerank, double *pagerank_old, int nb_nodes)
{
    /* update pagerank */
    for (int i = 0; i < nb_nodes; i++)
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
        /* if there is an edge from node to i */
        if (graph[node * nb_nodes + i])
        {
            sum++;
        }
    }

    return sum;
}