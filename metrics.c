/* Data Metrics routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include <string.h>

# include "global.h"

void get_fronts_from_file (int dtlz)
{
    FILE *gp_algorithm;
    FILE *gp_real_front;
    char str_algorithm[50],str_real_front[50];
    int i=0,j=0,k=0,p=3;
    char *token_algorithm,*token_real_front;
    const char s[2] = "\t";

    gp_algorithm=fopen("plot.out","rt");
    if (dtlz==1)
    	gp_real_front=fopen("real_front/DTLZ1-3-PF.txt","r");
    if (dtlz==2)
    	gp_real_front=fopen("real_front/DTLZ2-3-PF.txt","r");
    if (dtlz==3)
    	gp_real_front=fopen("real_front/DTLZ3-3-PF.txt","r");
    if (dtlz==4)
    	gp_real_front=fopen("real_front/DTLZ4-3-PF.txt","r");
    if (dtlz==5)
    	gp_real_front=fopen("real_front/DTLZ5-3-PF.txt","r");
    if (dtlz==6)
    	gp_real_front=fopen("real_front/DTLZ6-3-PF.txt","r");
    if (dtlz==7)
    	gp_real_front=fopen("real_front/DTLZ7-3-PF.txt","r");
    if( gp_real_front == NULL )
    {
      perror("Error while opening the real front file.\n");
      exit(EXIT_FAILURE);
    }
    if( gp_algorithm == NULL )
    {
      perror("Error while opening the algorithm front file.\n");
      exit(EXIT_FAILURE);
    }

    while (fgets(str_algorithm,50, gp_algorithm)!=NULL)
    {
	/*printf("%s\n",str_algorithm);*/
   	token_algorithm = strtok (str_algorithm,"\t");
		for (i=0;i<nobj;i++){	
		/*printf( " %f\n", atof(token_algorithm));*/
		igb_algorithm[i][j]=atof(token_algorithm);
		token_algorithm = strtok (NULL, "\t");
	   }
    j++;
    }
    while (fgets(str_real_front,50, gp_real_front)!=NULL)
    {
	/*printf("%s\n",str_real_front);*/
	token_real_front = strtok (str_real_front," ");
		for (i=0;i<nobj;i++){
		igb_real_front[i][k]=atof(token_real_front);
		/*printf( " %f\n", atof(token_real_front));*/
		token_real_front = strtok (NULL, " ");
	   }
    k++;
    }
   fclose(gp_algorithm);
   fclose(gp_real_front);
}
void IGD ()
{
	int i,j,k;
	double sum=0;
	double distance;
	double min_distance;
	double aux_distance;	
	
	get_maximum_value(1);
	get_minimum_value(1);

	get_normalized_front (0);
	get_normalized_front (1);

	
	for (k=0;k<popsize;k++)
	{
		min_distance=DBL_MAX;
		for (i=0;i<popsize;i++)
		{
			distance=0;
			for (j=0;j<nobj;j++)
			{
				distance+=pow(igb_real_front_normalized[j][k]-igb_algorithm_normalized[j][i],2);
			}
			aux_distance=sqrt(distance);
			if (aux_distance<min_distance)
			{
				min_distance=aux_distance;
			}
		}
	sum+=min_distance;
	}
	sum=sqrt(sum)/popsize;
	printf("the generational distance is %e\n",sum);
return;
}
void get_maximum_value(int front)
{
	int i,j;
	for (i=0;i<nobj;i++)
		maximum_value[i]=DBL_MIN;
	for (i=0;i<popsize;i++)
	{
		for (j=0;j<nobj;j++)
		{
			if (front)
			{	
				if (igb_real_front[j][i]>maximum_value[j])
				{
					maximum_value[j]=igb_real_front[j][i];
				}
			}
			else
			{
				if (igb_algorithm[j][i]>maximum_value[j])
				{
					maximum_value[j]=igb_algorithm[j][i];
				}
			}
		}
	}
	/*printf("The maximum values are: \n");
	for (i=0;i<nobj;i++)
		printf("%e\t",maximum_value[i]);
	printf("\n");*/
return;
}
void get_minimum_value(int front)
{
	int i,j;
	for (i=0;i<nobj;i++)
		minimum_value[i]=DBL_MAX;
	for (i=0;i<popsize;i++)
	{
		for (j=0;j<nobj;j++)
		{
			if (front)
			{
				if (igb_real_front[j][i]<minimum_value[j])
				{
					minimum_value[j]=igb_real_front[j][i];
				}
			}
			else
			{
				if (igb_algorithm[j][i]<minimum_value[j])
				{
					minimum_value[j]=igb_algorithm[j][i];
				}
			}
		}
	}
	/*printf("The minimum values are: \n");
	for (i=0;i<nobj;i++)
		printf("%e\t",minimum_value[i]);
	printf("\n");*/
return;
}
void get_normalized_front (int front)
{
	int i,j;
	for (i=0; i<popsize; i++)
	{
		for (j=0;j<nobj;j++)
		{
			if (front)
			{
				igb_real_front_normalized[j][i]=igb_real_front[j][i]-minimum_value[j]/(maximum_value[j]-minimum_value[j]);
			}
			else
			{
				igb_algorithm_normalized[j][i]=igb_algorithm[j][i]-minimum_value[j]/(maximum_value[j]-minimum_value[j]);
			}
		}
	}
return;
}
void hyper_volume()
{
	int i,j;
	
return;
}
