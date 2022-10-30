#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void read_ints(const char*, double**);
void update_hub(double*, double*, int);
void update_autority(double*, double*, int);

int main()
{
        double **graph;
        double *hub;
        double *autority;
        double *old_hub;
        double *old_autority;
         


        const int nb_nodes = 8;
        

        const int nb_iterations = 9;

        
        graph = calloc(nb_nodes, sizeof *graph);
        hub = calloc(nb_nodes, sizeof *hub);
        autority = calloc(nb_nodes, sizeof *autority);
        old_hub = calloc(nb_nodes, sizeof *old_hub);
        old_autority = calloc(nb_nodes, sizeof *old_autority);


        

        for (int i=0; i<nb_nodes; i++)
        {
                graph[i] = calloc(nb_nodes, sizeof *(graph[i]));
                autority[i]=1;
                hub[i]=1;
        }


        read_ints("test2.txt",graph);
        
        for(int k=0;k<nb_iterations;k++)
        {       
                
                update_hub(hub,old_hub,nb_nodes);
                update_autority(autority,old_autority,nb_nodes);
                
                printf("Iteration NB : %d\nHub :\t\tAutority :\n",k);
                for (int i=0; i<nb_nodes; i++)
                {
                        printf("%d -> %.5f , %d -> %.5f\n", i,old_hub[i],i,old_autority[i]);
                }
                printf("\n\n\n");	


                for (int i=0; i<nb_nodes; i++)
                {       
                        
                        for (int j=0; j<nb_nodes; j++)
                        {       
                                
                                if(graph[i][j])
                                {
                                        hub[i]+=old_autority[j];
                                }
                                if(graph[j][i])
                                {
                                        autority[i]+=old_hub[j];
                                }
                                
                        }
                        
                        
                }
                
                
        }

        free(graph);
        free(hub);
        free(autority);
        free(old_hub);
        free(old_autority);


        return 0;

}

void read_ints (const char* file_name, double** graph)
{
  FILE* file = fopen (file_name, "r");
  int i = 0;
  int j=0;

  
  while (!feof (file))
    {  
      
        fscanf (file, "%d", &i);
        fscanf (file, "%d", &j);      
        graph[i][j]=1;
    }
  fclose (file);        
}
void update_hub(double* hub, double* old_hub, int nb_nodes)
{       
        double sum=0;
        for(int i=0; i<nb_nodes; i++)
        {
                sum+=hub[i];
        }
        
        for (int i=0; i<nb_nodes; i++)
        {
                old_hub[i]=hub[i]/sum;
                hub[i]=0;
        }
}
void update_autority(double* autority, double* old_autority, int nb_nodes)
{       
        double sum=0;
        for(int i=0; i<nb_nodes; i++)
        {
                sum+=autority[i];
        }
        for (int i=0; i<nb_nodes; i++)
        {
                old_autority[i]=autority[i]/sum;
                autority[i]=0;
        }
}
