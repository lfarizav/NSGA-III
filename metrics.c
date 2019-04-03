/* Data Metrics routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include <string.h>

# include "global.h"

void invert_real_front(int dtlz)
{

    get_fronts_from_file (dtlz);
    FILE *real_front;
    real_front=fopen("real_front/IDTLZ1-3-PF.txt","w");
    int i,j,k;
    int n_var=5;
    int k1=n_var-nobj+1;
    int m;
    double g;
    for (i=0; i<popsize; i++)
    {
	for (j=0; j<nobj; j++)
	{
		g=0.0;
	    	for (m=0;m<nobj;m++)
	    	{
			g+=(ref_points[m][i]-0.5)*(ref_points[m][i]-0.5)-cos(20.0*PI*(ref_points[m][i]-0.5));	
	    	}
	    	g = 100*(k1+g);
		printf("g is %e\n",g);
		fprintf(real_front,"%e\t",(1.0+g)*0.5-igb_real_front[j][i]);
		
	}
	fprintf(real_front,"\n");
    }
    fflush(real_front);
    fclose(real_front);

}
void get_fronts_from_file (int dtlz)
{
    FILE *gp_algorithm;
    FILE *gp_real_front;
    char str_algorithm[1500],str_real_front[1500];
    int i=0,j=0,k=0,p=3;
    char *token_algorithm,*token_real_front;
    const char s[2] = "\t";
    gp_algorithm=fopen("plot.out","rt");
    if (dtlz==1)
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ1-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ1-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ1-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ1-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ1-15-PF.txt","r");		
    }
    else if (dtlz==2)
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ2-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ2-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ2-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ2-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ2-15-PF.txt","r");		
    }
    else if (dtlz==3)
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ3-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ3-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ3-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ3-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ3-15-PF.txt","r");		
    }
    else if (dtlz==4)
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ4-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ4-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ4-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ4-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ4-15-PF.txt","r");		
    }
    else if (dtlz==5)
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ5-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ5-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ5-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ5-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ5-15-PF.txt","r");		
    }
    else if (dtlz==7)
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ7-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ7-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ7-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ7-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ7-15-PF.txt","r");		
    }
    else if (dtlz==8)
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ1C-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ1C-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ1C-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ1C-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ1C-15-PF.txt","r");		
    }
    else
    {
	if (nobj==3)
    		gp_real_front=fopen("real_front/DTLZ6-3-PF.txt","r");
	if (nobj==5)
    		gp_real_front=fopen("real_front/DTLZ6-5-PF.txt","r");	
	if (nobj==8)
    		gp_real_front=fopen("real_front/DTLZ6-8-PF.txt","r");		
	if (nobj==10)
    		gp_real_front=fopen("real_front/DTLZ6-10-PF.txt","r");		
	if (nobj==15)	
    		gp_real_front=fopen("real_front/DTLZ6-15-PF.txt","r");	
    }

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

    while (fgets(str_algorithm,1500, gp_algorithm)!=NULL)
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
    while (fgets(str_real_front,1500, gp_real_front)!=NULL)
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
   IGDfrontsize=k;
}
double IGD (population *pop)
{
	int i,j,k;
	double sum=0;
	double distance;
	double min_distance;
	double aux_distance;	

	/*get_maximum_value(1);
	get_minimum_value(1);

	get_normalized_front (0);
	get_normalized_front (1);*/
        if (dtlz<16)
        	get_fronts_from_file (dtlz);

	for (k=0;k<factorial+factorial_inside;k++)
	{
		min_distance=DBL_MAX;
		for (i=0;i<factorial+factorial_inside;i++)
		{
			distance=0;
			for (j=0;j<nobj;j++)
			{

				/*printf("[%d] Real front %e, NSGA-III front %e\n",i,igb_real_front[j][i],igb_algorithm[j][i]);*/
				if (dtlz<16)
					distance+=pow(igb_real_front[j][k]-igb_algorithm[j][i],2);
				else
				{
					distance+=pow(DTLZ[j][k]-pop->ind[i].obj[j],2);
					/*printf("real %e, algorithm %e\n",DTLZ[j][i],pop->ind[i].obj[j]);*/
				}
					
			}
			aux_distance=sqrt(distance);
			if (aux_distance<min_distance)
			{
				min_distance=aux_distance;
			}
		}
		sum+=min_distance;
	}
	/*printf("IGD %e\n",sum/(factorial+factorial_inside));*/
return sum/(factorial+factorial_inside);
}
void get_maximum_value(int front)
{
	int i,j;
	for (i=0;i<nobj;i++)
		maximum_value[i]=DBL_MIN;
	for (i=0;i<IGDfrontsize;i++)
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
	for (i=0;i<IGDfrontsize;i++)
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
	for (i=0; i<IGDfrontsize; i++)
	{
		for (j=0;j<nobj;j++)
		{
			/*printf("real %e, algorithm %e\n",igb_real_front[j][i],igb_algorithm[j][i]);*/
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
