/*
 *
 * @author : Cybermade
 * gcc -Wall -o pagerank_algorithm pagerank_algorithm.c
 * ./pagerank_algorithm <input_file> <nb_nodes> <nb_iterations> 
 *
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int number_of_neighbors(int **, int , int);
void update_pagerank(double*, double*, int);
void read_ints(const char *file_name, int **graph)
{
    FILE *file = fopen(file_name, "r");
    int i = 0;
    int j = 0;

    while (!feof(file))
    {

        fscanf(file, "%d", &i);
        fscanf(file, "%d", &j);
        graph[i][j] = 1;
    }
    fclose(file);
}
int main(int argc, char *argv[])
{
    int **graph;
    char *nbnodes = argv[2];
    int nb_nodes = atoi(nbnodes); /* number of nodes */
    double dumping_factor = 0.85;
    double *pagerank;
    double *pagerank_old;
    char *nbiterations = argv[3];
    int nb_iterations = atoi(nbiterations);  /* number of iterations to run the algorithm */

    graph = calloc(nb_nodes, sizeof *graph);
    pagerank = calloc(nb_nodes, sizeof (double));
    pagerank_old = calloc(nb_nodes, sizeof (double));
    
    for (int i = 0; i < nb_nodes; i++)
    {
        graph[i] = calloc(nb_nodes, sizeof *(graph[i]));
        
    }
    read_ints(argv[1], graph);
    for (int i = 0; i < nb_nodes; i++)
    {
        for (int j = 0; j < nb_nodes; j++)
        {
            printf("%d ", graph[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < nb_nodes; i++)
    {
        pagerank[i] = 1/(double)nb_nodes;
    }
    
    for (int k =0;k<=nb_iterations;k++)
    {   printf("Iteration %d\n", k);
        for (int i = 0; i < nb_nodes; i++)
        {
        printf("%d -> %.5f\n", i, pagerank[i]);
        }
        double sum = 0;
        for (int i = 0; i < nb_nodes; i++)
        {
            sum+= pagerank[i];
            
        }
        printf("Sum : %.5f (Should always be 1)\n\n", sum);
        update_pagerank(pagerank, pagerank_old, nb_nodes);
    
    for(int i=0;i<nb_nodes;i++)
    {   pagerank[i] = (1 - dumping_factor) /nb_nodes;
        for(int j=0;j<nb_nodes;j++)
        {
            if(graph[j][i])
            {   
                pagerank[i] += (dumping_factor)*(pagerank_old[j]/(double)number_of_neighbors(graph, j, nb_nodes));
            }
        }
        
        
    }
    
    }
    
    for (int i = 0; i < nb_nodes; i++)
    {
        free(graph[i]);
    }
    free(graph);
    free(pagerank);
    return 0;
}
void update_pagerank(double *pagerank, double *pagerank_old, int nb_nodes)
{
    for(int i = 0; i < nb_nodes; i++)
    {
        pagerank_old[i] = pagerank[i];
        pagerank[i] = 0;
    }
}
int number_of_neighbors(int **graph, int node, int nb_nodes)
{
    int sum = 0;
    for (int i = 0; i < nb_nodes; i++)
    {   
        if(graph[node][i])
        {
            sum++;
        }
        
    }
    return sum;
}