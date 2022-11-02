#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
int main()
{

    int **graph;
    double *hub;
    double *autority;

    double norm = 0;

    const int nb_nodes = 8;

    const int nb_iterations = 10;

    graph = calloc(nb_nodes, sizeof *graph);
    hub = calloc(nb_nodes, sizeof (double));
    autority = calloc(nb_nodes, sizeof (double));
    for (int i = 0; i < nb_nodes; i++)
    {
        graph[i] = calloc(nb_nodes, sizeof *(graph[i]));
        autority[i] = 1;
        hub[i] = 1;
    }
    read_ints("test2.txt", graph);

    printf("Iteration NB : 0\nHub :\t\tAutority :\n");
        for (int i = 0; i < nb_nodes; i++)
        {
            printf("%d -> %.5f , %d -> %.5f\n", i, hub[i], i, autority[i]);
        }
    for (int k = 1; k <= nb_iterations; k++)
    {

        norm = 0;

        for (int i = 0; i < nb_nodes; i++)
        {
            autority[i] = 0;
            for (int j = 0; j < nb_nodes; j++)
            {

                if (graph[j][i])
                {
                    autority[i] += hub[j];
                }
            }
            norm += pow(autority[i], 2);
        }
        norm = sqrt(norm);
        for (int i = 0; i < nb_nodes; i++)
        {
            autority[i] /= norm;
        }


        
        norm = 0;
        for (int i = 0; i < nb_nodes; i++)
        {
            hub[i] = 0;
            for (int j = 0; j < nb_nodes; j++)
            {

                if (graph[i][j])
                {
                    hub[i] += autority[j];
                }
            }
            norm += pow(hub[i], 2);
        }
        norm = sqrt(norm);
        for (int i = 0; i < nb_nodes; i++)
        {
            hub[i] /= norm;
        }
        printf("Iteration NB : %d\nHub :\t\tAutority :\n", k);
        for (int i = 0; i < nb_nodes; i++)
        {
            printf("%d -> %.5f , %d -> %.5f\n", i, hub[i], i, autority[i]);
        }
        printf("\n\n\n");
    }

    /*
    G := set of pages
for each page p in G do
    p.auth = 1 // p.auth is the authority score of the page p
    p.hub = 1 // p.hub is the hub score of the page p
for step from 1 to k do // run the algorithm for k steps
    norm = 0
    for each page p in G do  // update all authority values first
        p.auth = 0
        for each page q in p.incomingNeighbors do // p.incomingNeighbors is the set of pages that link to p
            p.auth += q.hub
        norm += square(p.auth) // calculate the sum of the squared auth values to normalise
    norm = sqrt(norm)
    for each page p in G do  // update the auth scores
        p.auth = p.auth / norm  // normalise the auth values
    norm = 0
    for each page p in G do  // then update all hub values
        p.hub = 0
        for each page r in p.outgoingNeighbors do // p.outgoingNeighbors is the set of pages that p links to
            p.hub += r.auth
        norm += square(p.hub) // calculate the sum of the squared hub values to normalise
    norm = sqrt(norm)
    for each page p in G do  // then update all hub values
        p.hub = p.hub / norm   // normalise the hub values*/

    return 0;
}