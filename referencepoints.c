/* Data initializtion routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include <string.h>

# include "global.h"

/* Function to generate reference points*/
/*step = nobj,1/(double)(numberpointperdim-1),factorial*/
int generate_DTLZ1 (int nobj_for, double step)
{
    int i,count;
    count=recursive_for_DTLZ1 (nobj_for,step,0,0);
    /*printf("There are %d reference points \n\n",count);*/
    return count;
}
int recursive_for_DTLZ1 (int nobj_for, double step, int count, int i)
{
	int j;
	
        if (count==factorial )
		return 1;
    	if (count>factorial)
	{
		printf("Error, count > factorial %d, then exit!, check recursive_for function\n",factorial);
		exit(-1);
	}
	for (j=0;j<numberpointperdim-i && count<factorial;j++)
	{	
		int nobj_for_temp=nobj_for;
		index[nobj-1-nobj_for_temp]=j;

		if (nobj_for_temp==1)
		{
			int temp=0;
			int k;

			int m;
			int n_var=5;
		    	int k1=n_var-nobj+1;
		    	double g=0.0;
		    	for (m=n_var-k1;m<n_var;m++)
		    	{	
				g+=(index[m]-0.5)*(index[m]-0.5)-cos(20.0*PI*(index[m]-0.5));	
		    	}
		    	g = 100*(k1+g);
			for (k=0;k<nobj-1;k++)
			{
				DTLZ[k][count]=0.5-0.5*step*index[k];
				temp+=index[k];
			}
			DTLZ[k][count]=0.5-0.5*(1-(step*temp));
			count++;
		}
		if (nobj_for_temp>1 && count<factorial)
		{
			nobj_for_temp--;
			count=recursive_for_DTLZ1(nobj_for_temp,step,count,i+j);
		}
	}
	
	return count;
}

int generate_ref_points (int nobj_for, double step)
{
    int i,count;
    	count=recursive_for (nobj_for,step,0,0);
	/*printf("There are %d reference points \n\n",count);*/
    return count;
}
int recursive_for (int nobj_for, double step, int count, int i)
{
	int j;
	
        if (count==factorial )
		return 1;
    	if (count>factorial)
	{
		printf("Error, count > factorial %d, then exit!, check recursive_for function\n",factorial);
		exit(-1);
	}
	for (j=0;j<numberpointperdim-i && count<factorial;j++)
	{	
		int nobj_for_temp=nobj_for;
		index[nobj-1-nobj_for_temp]=j;

		if (nobj_for_temp==1)
		{
			int temp=0;
			int k;
			for (k=0;k<nobj-1;k++)
			{
				ref_points[k][count]=step*index[k];
				temp+=index[k];
			}
			ref_points[k][count]=1-step*(temp);
			count++;
		}
		if (nobj_for_temp>1 && count<factorial)
		{
			nobj_for_temp--;
			count=recursive_for(nobj_for_temp,step,count,i+j);
		}
	}
	
	return count;
}
/* Function to generate reference points in the inside layer (two layer case)*/
int generate_ref_points_inside (int nobj_for, double step)/*nobj,1/(double)(numberpointperdim-1),factorial*/
{
    int count;

    count=recursive_for_inside (nobj_for,step,0,0);
    /*printf("There are %d inside reference points \n\n",count);*/
    return count;
}
int recursive_for_inside (int nobj_for, double step, int count, int i)
{
	int j;
        if (count==factorial_inside )
		return 1;
    	if (count>factorial_inside)
	{
		printf("error, count > factorial_inside %d, then exit!, check recursive_for function\n",factorial_inside);
		exit(-1);
	}
	for (j=0;j<numberpointperdim_inside-i && count<factorial_inside;j++)
	{	
		int nobj_for_temp=nobj_for;
		index[nobj-1-nobj_for_temp]=j;

		if (nobj_for_temp==1)
		{
			int temp=0;
			int k;
			for (k=0;k<nobj-1;k++)
			{
				ref_points[k][count+factorial]=(1-step)*step*index[k]+(double)step/nobj;
				temp+=index[k];
			}
			ref_points[k][count+factorial]=(1-step)*(1-step*(temp))+(double)step/nobj;
			count++;
		}
		if (nobj_for_temp>1 && count<factorial_inside)
		{
			nobj_for_temp--;
			count=recursive_for_inside(nobj_for_temp,step,count,i+j);
		}
	}
	
	return count;
}
/* Function to generate adaptive reference points*/
int generate_adaptive_ref_points (int nobj_for, double step)/*nobj,1/(double)(numberpointperdim-1),factorial*/
{
    int i;
    printf("Adaptive vectors are:\n");
    return recursive_for_adaptive (nobj_for,step,0,0);
}
int recursive_for_adaptive (int nobj_for, double step, int count, int i)
{
	int j;
	
        if (count==factorial_adaptive )
		return 1;
    	if (count>factorial_adaptive)
	{
		printf("error, count > factorial_adaptive = %d, then exit!, check recursive_for function\n", factorial_adaptive);
		exit(-1);
	}
	for (j=0;j<2-i && count<factorial_adaptive;j++)
	{	
		int nobj_for_temp=nobj_for;
		index[nobj-1-nobj_for_temp]=j;

		if (nobj_for_temp==1)
		{
			int temp=0;
			int k;
			for (k=0;k<nobj-1;k++)
			{
				minimum_amount_ref_points[k][count]=(1/(double)(numberpointperdim-1))*step*index[k];
				printf("%e\t",minimum_amount_ref_points[k][count]);					
				temp+=index[k];
			}
			
			minimum_amount_ref_points[k][count]=(1/(double)(numberpointperdim-1))*(1-step*(temp));
			printf("%e\n",minimum_amount_ref_points[k][count]);
			count++;
		}
		if (nobj_for_temp>1 && count<factorial_adaptive)
		{
			nobj_for_temp--;
			count=recursive_for_adaptive(nobj_for_temp,step,count,i+j);
		}
	}
	
	return count;
}
long fact (int x)
{
     long int f=1;
     int i;
     for(i=1;i<=x;i++){
       f=f*i;
     }
     return(f);
}
void square_refpoints()
{
	int i,j;
	for (i=0;i<factorial+factorial_inside;i++)
	{
		for (j=0;j<nobj;j++)
		{
			ref_points[j][i]=sqrt(ref_points[j][i]);
		}
	}
return;
}
void refpoints_normalized()
{
	int i,j;
	min_refpoints();
	max_refpoints();
	for (i=0;i<factorial+factorial_inside;i++)
	{
		for (j=0;j<nobj;j++)
		{
			ref_points_normalized[j][i]=(ref_points[j][i]-min_ref_points[j])/(max_ref_points[j]-min_ref_points[j]);
		}
	}
return;
}
void min_refpoints()
{
	int i,j;
	for (j=0;j<nobj;j++)
	{
		min_ref_points[j]=DBL_MAX;
	}
	for (i=0;i<factorial+factorial_inside;i++)
	{
		for (j=0;j<nobj;j++)
		{
			if (ref_points[j][i]<min_ref_points[j])
				min_ref_points[j]=ref_points[j][i];
		}
	}
	printf("The minimun values of reference points are:\n");
	for (j=0;j<nobj;j++)
	{
		printf("%e\n",min_ref_points[j]);
	}
	printf("\n");
return;
}
void max_refpoints()
{
	int i,j;
	for (j=0;j<nobj;j++)
	{
		min_ref_points[j]=0;
	}
	for (i=0;i<factorial+factorial_inside;i++)
	{
		for (j=0;j<nobj;j++)
		{
			if (ref_points[j][i]>max_ref_points[j])
				max_ref_points[j]=ref_points[j][i];
		}
	}
	printf("The maximum values of reference points are:\n");
	for (j=0;j<nobj;j++)
	{
		printf("%e\n",max_ref_points[j]);
	}
	printf("\n");
return;
}
int get_supplied_reference_points_from_file (char *filename)
{
    FILE *supplied_ref_points;
    char str_supplied_ref_points[50];
    int i;
    int j=0;
    char *token_supplied_ref_points;
    double s;

    supplied_ref_points=fopen(filename,"rt");

    if( supplied_ref_points == NULL )
    {
      perror("Error while opening supplied reference points file.\n");
      exit(EXIT_FAILURE);
    }

    while (fgets(str_supplied_ref_points,50, supplied_ref_points)!=NULL)
    {
	/*printf("%s\n",str_supplied_ref_points);*/
   	token_supplied_ref_points = strtok (str_supplied_ref_points,"\t");
		for (i=0;i<nobj;i++)
		{			
			s=atof(token_supplied_ref_points);
			printf( "%f\n",s);
			token_supplied_ref_points = strtok (NULL, "\t");
	    	}
    j++;
    }
   fclose(supplied_ref_points);
return j;
}

int create_adaptive_refpoints()
{
	/*This function take information of the updated niche count (Next size population size(P(t+1))=N )*/
	int i,j,k,l,m,n;
	int crowded_rho=1;
	int temp=0;
	int temp_j=0;
	int temp_nobj;
	int temp_nobj_adaptive;
	int temp_outside_first_quadrant;
	int temp_adaptive_refpoint_repeated;

	/*Organize original Das and Dennis refpoints from highest to lowest rho_St*/
    	sort_all_refpoints_by_rho_index();
	/*find_useless_refpoints_index finds index for useless_refpoints and initialize usefull_refpoint_number*/
        find_useless_usefull_refpoint_index();

	printf("\nFinding the number of useless and usefull reference points\n_______________________________________________\n");
	printf("The number of useless reference points is %d, the number of usefull reference points is %d, Total %d\n", useless_refpoint_number, usefull_refpoint_number,useless_refpoint_number+usefull_refpoint_number);
	/*Check rho_St>1 just for usefull refpoints (Original Das an Dennis reference points)*/
	for (i=0;i<usefull_refpoint_number;i++)
	{
		printf("i %d, sort_all_refpoint_index[%d] %d, rho_St[%d] %d, rho_Fl[%d] %d\n",i,i,sort_all_refpoint_index[i],sort_all_refpoint_index[i],rho_St[sort_all_refpoint_index[i]],sort_all_refpoint_index[i],rho_Fl[sort_all_refpoint_index[i]]);
		for (j=0;j<factorial_adaptive;j++)
		{
		    /*That is true because refpoints are sorted from the highest to the lowest value. So, 
	    	    it requires to do only usefull_refpoints_number times*/
		    if (rho_St[sort_all_refpoint_index[i]]>crowded_rho)
		    {
		      if (adaptive_nsga==1)
		      {
			printf("Adaptive refpoints:\n");
			for (k=0;k<nobj;k++)
			{
				/*Generation and Translation of new reference points to the crowded reference point*/
				adaptive_refpoints[k][temp_j]=minimum_amount_ref_points[k][j]+ref_points[k][sort_all_refpoint_index[i]]-(1/(double)(numberpointperdim-1))/nobj;
				printf("%e\t",adaptive_refpoints[k][temp_j]);
			}
			printf("\n");
			/*Verification:reference points adaptive neither outside the first quadrant nor repeated*/
			for (l=0;l<factorial+factorial_inside+last_gen_adaptive_refpoints_number;l++)
			{
				temp_nobj=0;
				temp_outside_first_quadrant=0;
				for (k=0;k<nobj;k++)
				{
					if (fabs(adaptive_refpoints[k][temp_j]-ref_points[k][l])<1e-15)
					{
						temp_nobj++;	
					}
					if (adaptive_refpoints[k][temp_j]<0)
					{
						temp_outside_first_quadrant=1;
						printf("adaptive_refpoints[%d][%d]=%e,outside the first quadrant\n",k,temp_j,adaptive_refpoints[k][temp_j]);
					}

				}
				if (temp_nobj==nobj || temp_outside_first_quadrant)
				{
					printf("adaptive_refpoint [%d] repeated? %d, temp_outside_first_quadrant? %d\n",
					temp_j,temp_nobj==nobj,temp_outside_first_quadrant);
					break;/*l=factorial+factorial_inside+last_gen_adaptive_refpoints_number;*/
				}
			}
			/*Be sure there is not repeated the new adaptive reference point, before add it*/	
			for (m=0;m<temp_j;m++)
			{
				temp_nobj_adaptive=0;
				temp_adaptive_refpoint_repeated=0;
				for (k=0;k<nobj;k++)
				{
					if (fabs(adaptive_refpoints[k][temp_j]-adaptive_refpoints[k][m])<1e-15)
						temp_nobj_adaptive++;
				}
				if (temp_nobj_adaptive==nobj)
				{
					temp_adaptive_refpoint_repeated=1;
					printf("Adaptive reference point repeated\n");
					break;
				}
			}	
			if (!(temp_nobj==nobj || temp_outside_first_quadrant || temp_adaptive_refpoint_repeated))
			{
				printf("adaptive_refpoint [%d] added, after doing 3 checks\n",temp_j);
				temp_j++;
			}
		      }
		      if (adaptive_nsga==2)
		      {
			double translation[nobj][nobj];
			for (l=0;l<nobj;l++)
			{
				
				for (k=0;k<nobj;k++)
				{
					translation[k][l]=0.5*minimum_amount_ref_points[k][l]-(1/(double)(numberpointperdim-1))/(2*nobj);
				}
			}
			printf("Adaptive points translated to the crowded reference point\n");
			for (k=0;k<nobj;k++)
			{
				for (l=0;l<nobj;l++)
				{
					adaptive_refpoints[l][temp_j]=translation[l][k]-translation[l][j]+ref_points[l][sort_all_refpoint_index[i]];
				}
				/*Verification:reference points adaptive neither outside the first quadrant nor repeated*/
				/*It either repeated of outside the first quadrant, assign the vector 0*/
				for (l=0;l<factorial+factorial_inside+last_gen_adaptive_refpoints_number;l++)
				{
					temp_nobj=0;
					temp_outside_first_quadrant=0;
					for (m=0;m<nobj;m++)
					{
						if (fabs(adaptive_refpoints[m][temp_j]-ref_points[m][l])<1e-15)
						{
							temp_nobj++;	
						}
						if (adaptive_refpoints[m][temp_j]<0)
						{
							temp_outside_first_quadrant=1;
							printf("adaptive_refpoints[%d][%d]=%e,outside the first quadrant\n",m,temp_j,adaptive_refpoints[m][temp_j]);
						}

					}
					if (temp_nobj==nobj || temp_outside_first_quadrant)
					{
						printf("adaptive_refpoint [%d] repeated? %d, temp_outside_first_quadrant? %d\n",
						temp_j,temp_nobj==nobj,temp_outside_first_quadrant);
						l=factorial+factorial_inside+last_gen_adaptive_refpoints_number;
					}
				}
				/*Be sure there is not repeated adaptive refpoints, before add*/
				for (m=0;m<temp_j;m++)
				{
					temp_nobj_adaptive=0;
					temp_adaptive_refpoint_repeated=0;
					for (n=0;n<nobj;n++)
					{
						if (fabs(adaptive_refpoints[n][temp_j]-adaptive_refpoints[n][m])<1e-15)
							temp_nobj_adaptive++;
					}
					if (temp_nobj_adaptive==nobj)
					{
						temp_adaptive_refpoint_repeated=1;
						printf("Adaptive reference point repeated\n");
						break;
					}
				}	
				if (!(temp_nobj==nobj || temp_outside_first_quadrant || temp_adaptive_refpoint_repeated))
				{
					printf("adaptive_refpoint [%d] added, after doing 3 checks\n",temp_j);
					temp_j++;
				}
			}
		      }
		    }
		}
	}
	printf("The number of adaptive refpoints neither repeated nor outside the first quadrant is %d\n",temp_j);
	printf("Vacancies to fill, instead of original reference points with (rho_St = 0) = %d\n",useless_refpoint_number);

return temp_j;
}

void add_adaptive_refpoints_to_ref_points()
{
	int i,j,k,l,m;
	int temp;
	int temp_nobj;
        /*create new adaptive reference points */
	/*Here adaptive_refpoint_number is returned*/
	printf("After niching, We add new adaptive reference points\n");
   	adaptive_refpoint_number=create_adaptive_refpoints();
        /*replace useless refpoints with adaptive refpoints*/
	printf("Following values of rho_St and rho_Fl must be 0\n");
	printf("adaptive refpoints number  %d\n",adaptive_refpoint_number);
	for(i=0;i<adaptive_refpoint_number;i++)
	{
		printf("%d, rho_St[%d] %d,rho_Fl[%d] %d,rho[%d] %d\n",i,useless_refpoint_index[i],rho_St[useless_refpoint_index[i]],useless_refpoint_index[i],rho_Fl[useless_refpoint_index[i]],useless_refpoint_index[i],rho[useless_refpoint_index[i]]);
	}
	printf("factorial=%d, factorial_inside=%d, last_gen_adaptive_refpoints_number=%d, adaptive_refpoint_number=%d\n",factorial,factorial_inside,last_gen_adaptive_refpoints_number,adaptive_refpoint_number);
	printf("The adaptive refpoints are:\n");
	for(i=0;i<adaptive_refpoint_number;i++)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			printf("%e\t",adaptive_refpoints[k][i]);
		}
		printf("\n");
	}
	printf("The reference points further Das and Dennis are:\n");
	/*Be sure there is not repeated the new adaptive reference point, before add it*/
	for(j=factorial+factorial_inside+last_gen_adaptive_refpoints_number;j<factorial+factorial_inside+last_gen_adaptive_refpoints_number+adaptive_refpoint_number;j++)
	{
		printf("%d\t",j);
		for (k=0;k<nobj;k++)
		{
			ref_points[k][j]=0.0;
			printf("%e\t",ref_points[k][j]);
		}
		printf("\n");
	}
	printf("\nThe adaptive refpoints are:\n");
	for(i=0;i<adaptive_refpoint_number;i++)
	{
		printf("%d\t",i+factorial+factorial_inside+last_gen_adaptive_refpoints_number);
		for (k=0;k<nobj;k++)
		{
			ref_points[k][factorial+factorial_inside+last_gen_adaptive_refpoints_number+i]=adaptive_refpoints[k][i];
			printf("%e\t",adaptive_refpoints[k][i]);
		}
		printf("\n");
	}
        display_refpoints ();
return;
}
int delete_adaptive_refpoints(int archieve_size, int front_size, population *selection_pop, population *new_pop, int generation)
{	
	/*If there are not adaptive reference points added, then go directly to the next generation*/
	if (adaptive_refpoint_number==0)
		return 0;
	int i,j,k,l;
	int temp_adaptive_number;
	int temp_nobj;
	int is_adaptive_refpoint_already_included_in_das_and_dennis_refpoints=0;

	printf("factorial %d, adaptive_ref_points_inserted %d, adaptive_ref_points_inserte_per_generation %d, factorial_inside %d, last_gen_adaptive_refpoints_number %d, adaptive_refpoint_number %d\n",factorial,adaptive_ref_points_inserted,adaptive_ref_points_inserted_per_generation,factorial_inside, last_gen_adaptive_refpoints_number, adaptive_refpoint_number);
    	/* Routine to find min from functions*/
    	for (j=0; j<popsize; j++)
    	{
      	  find_min_from_functions(&(new_pop->ind[j]),j,1);/*ok*/
	  find_max_from_functions(&(new_pop->ind[j]),j,1);/*ok*/
        }
    	/*Routine that substracts the zmin to solutions*/
    	for (k=0; k<popsize; k++)
    	{
        	obj_minus_zmin(&(new_pop->ind[k]));/*ok*/
    	}
    	find_extreme_points(new_pop, popsize);/*only consider the individuals in the first front*/

    	/*normalize solutions*/
    	for (l=0; l<popsize; l++)
    	{
        	normalized_objective_function(&(new_pop->ind[l]));/*ok*/
    	}
	/*If there are adaptive reference points, then associate them to the population*/
    	for (l=0; l<factorial+factorial_inside+last_gen_adaptive_refpoints_number+adaptive_refpoint_number; l++)
    	{
		rho_St[l]=0;
		rho_Fl[l]=0;
		rho[l]=0;
	}
    	for (l=0; l<popsize; l++)
    	{
		associate(&new_pop->ind[l],&new_pop->ind[l],l,10*popsize,0,factorial+factorial_inside+last_gen_adaptive_refpoints_number+adaptive_refpoint_number);
	}
	printf("Printing rho after adding adaptive refpoints\n");
    	int temp_rho_St_total=0;
    	int temp_rho_Fl_total=0;
	for (l=0; l<factorial+factorial_inside+last_gen_adaptive_refpoints_number+adaptive_refpoint_number; l++)
	{
		temp_rho_St_total+=rho_St[l];
		temp_rho_Fl_total+=rho_Fl[l];
		printf("rho_St[%d] %d, rho_Fl[%d] %d\t",l,rho_St[l],l,rho_Fl[l]);
		for (k=0;k<nobj;k++)
		{
			printf("%e\t",ref_points[k][l]);
		}
		printf("\n");
	}
    	printf("After adding adaptive refpoints: temp_rho_St_total %d, temp_rho_Fl_total %d\n",temp_rho_St_total,temp_rho_Fl_total);
	/*Add to the next generation all initial reference points plus adaptive reference points with rho_St=1*/
	check_adaptive_refpoints_inclusion_number(generation);
	printf("factorial %d, factorial_inside %d, last_gen_adaptive_refpoints_number %d, adaptive_ref_points_inserted_per_generation %d\n",factorial,factorial_inside,last_gen_adaptive_refpoints_number,adaptive_ref_points_inserted_per_generation);
	temp_adaptive_number=0;
	/*for (i=factorial-adaptive_ref_points_inserted_per_generation+factorial_inside+last_gen_adaptive_refpoints_number;i<factorial-adaptive_ref_points_inserted_per_generation+factorial_inside+last_gen_adaptive_refpoints_number+adaptive_refpoint_number;i++)
	{
		if (rho_St[i]==1)
		{
			printf("%d\t",i);
			for (k=0;k<nobj;k++)
			{
				printf("%e\t",adaptive_refpoints[k][i-(factorial+factorial_inside+last_gen_adaptive_refpoints_number-adaptive_ref_points_inserted_per_generation)]);
			}
			printf("\n");
			for (j=0;j<adaptive_ref_points_inserted_per_generation;j++)
			{
				temp_nobj=0;
				is_adaptive_refpoint_already_included_in_das_and_dennis_refpoints=0;
				printf("Checking refpoints already in DasandDennis refponits\n");
				for (k=0;k<nobj;k++)
				{
					printf("(%e,%e)\n",ref_points[k][factorial-adaptive_ref_points_inserted_per_generation+factorial_inside+j],adaptive_refpoints[k][i-(factorial+factorial_inside+last_gen_adaptive_refpoints_number-adaptive_ref_points_inserted_per_generation)]);
					if(fabs(ref_points[k][factorial-adaptive_ref_points_inserted_per_generation+factorial_inside+j]-adaptive_refpoints[k][i-(factorial+factorial_inside+last_gen_adaptive_refpoints_number-adaptive_ref_points_inserted_per_generation)])<1e-15)
						temp_nobj++;	
				}
				if (temp_nobj==nobj)
				{
					is_adaptive_refpoint_already_included_in_das_and_dennis_refpoints=1;
					break;
				}
			}
			printf("i %d, rho_St[%d] %d\n",i,i,rho_St[i]);
			printf("%d\t",i);
			if (!is_adaptive_refpoint_already_included_in_das_and_dennis_refpoints)
			{
				for (k=0;k<nobj;k++)
				{
					ref_points[k][factorial+factorial_inside+temp_adaptive_number]=adaptive_refpoints[k][i-(factorial+factorial_inside+last_gen_adaptive_refpoints_number-adaptive_ref_points_inserted_per_generation)];	
					printf("%e\t",ref_points[k][factorial+factorial_inside+temp_adaptive_number]);
				}
				printf("\n");
				temp_adaptive_number++;
			}
		}
		
	}*/
	printf("is_adaptive_refpoint_already_included_in_das_and_dennis_refpoints=%d\n",is_adaptive_refpoint_already_included_in_das_and_dennis_refpoints);
	printf("Visualization of last generation adaptive refpoints, before next generation. Adding %d adaptive refpoints\n", temp_adaptive_number);
	for (i=factorial+factorial_inside;i<factorial+factorial_inside+temp_adaptive_number;i++)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			printf("%e\t",ref_points[k][i]);
		}
		printf("\n");
	}
	adaptive_refpoint_number=0;
return temp_adaptive_number;
}

void check_adaptive_refpoints_inclusion_number(int generation)
{
    int i,j,k;
    int temp_nobj;
    int temp_j=0;
    adaptive_ref_points_inserted_per_generation=0;
    printf("\n");
    for (i=factorial+factorial_inside+last_gen_adaptive_refpoints_number;i<factorial+factorial_inside+last_gen_adaptive_refpoints_number+adaptive_refpoint_number;i++)
    {
	if (rho_St[i]>0)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			printf("%e\t",ref_points[k][i]);
		}
		printf("\n");
		for (j=0;j<elegible_adaptive_ref_points_to_be_fixed_number;j++)
		{
			temp_nobj=0;
			printf("(adaptive_ref_point_settled[k][%d],ref_point[k][%d])\n",j,i);
			for (k=0;k<nobj;k++)
			{
				if (fabs(adaptive_ref_points_settled[k][j]-ref_points[k][i])<1e-15)
				{
					temp_nobj++;
				}
				printf("(%e-%e)=%e\t",adaptive_ref_points_settled[k][j],ref_points[k][i],fabs(adaptive_ref_points_settled[k][j]-ref_points[k][i]));
			}
			printf("\n");
			if (temp_nobj==nobj)
			{
				printf("%d, the same!!\n",j);
				if (rho_St[i]!=last_rho_St[j])
				{
					last_rho_St[j]=rho_St[i];
					last_generation_associated_rho_St[j]=generation;
					adaptive_ref_points_settled_number[j]=0;
					printf("rho_St %d,last_rho_St %d \n", rho_St[i],last_rho_St[j]);
				}
				temp_j=j;
				break;
			}
		}
		if (temp_nobj==nobj)
		{
			printf("Settled adaptive reference point is already included, increase the ocurrence number\n");
			adaptive_ref_points_settled_number[temp_j]=adaptive_ref_points_settled_number[temp_j]+1;
			printf("elegible_adaptive_ref_points_to_be_fixed_number %d\n",elegible_adaptive_ref_points_to_be_fixed_number);
			if (adaptive_ref_points_settled_number[temp_j]>=3 &&(generation-last_generation_associated_rho_St[temp_j])>=10)
			{
				printf("New reference point added to Initial Das and Dennis reference points\n");
				for (k=0;k<nobj;k++)
				{
					ref_points[k][useless_refpoint_index[temp_i]]=adaptive_ref_points_settled[k][temp_j];
					adaptive_ref_points_settled[k][temp_j]=adaptive_ref_points_settled[k][elegible_adaptive_ref_points_to_be_fixed_number-1];
					printf("%e\t",ref_points[k][factorial+factorial_inside]);
				
				}
				printf("\n");
				adaptive_ref_points_settled_number[temp_j]=adaptive_ref_points_settled_number[elegible_adaptive_ref_points_to_be_fixed_number-1];
				last_rho_St[temp_j]=rho_St[elegible_adaptive_ref_points_to_be_fixed_number-1];
				adaptive_ref_points_settled_number[elegible_adaptive_ref_points_to_be_fixed_number-1]=0;
				elegible_adaptive_ref_points_to_be_fixed_number--;
				factorial++;
				adaptive_ref_points_inserted++;
				adaptive_ref_points_inserted_per_generation++;
			}

		}
		else
		{
			last_rho_St[elegible_adaptive_ref_points_to_be_fixed_number]=rho_St[i];
			printf("New settled adaptive reference point is included, the ocurrence number is \t");
			for (k=0;k<nobj;k++)
			{
				adaptive_ref_points_settled[k][elegible_adaptive_ref_points_to_be_fixed_number]=ref_points[k][i];
			}
			adaptive_ref_points_settled_number[elegible_adaptive_ref_points_to_be_fixed_number]=adaptive_ref_points_settled_number[elegible_adaptive_ref_points_to_be_fixed_number]+1;
			last_generation_associated_rho_St[elegible_adaptive_ref_points_to_be_fixed_number]=generation;
			printf("%d\n",adaptive_ref_points_settled_number[elegible_adaptive_ref_points_to_be_fixed_number]);
			elegible_adaptive_ref_points_to_be_fixed_number++;
		}
	}
    }
    printf("Printing the adaptive_ref_points_settled\n");
    for (i=0;i<elegible_adaptive_ref_points_to_be_fixed_number;i++)
    {
	printf("i %d, ocurrence %d, last_rho_St %d, last_generation_associated_rho_St %d, generation %d\t",i,adaptive_ref_points_settled_number[i],last_rho_St[i],last_generation_associated_rho_St[i],generation);
	for (j=0;j<nobj;j++)
	{
		printf("\t %e",adaptive_ref_points_settled[j][i]);
	}
	printf("\n");
    }
return;
}
/*This function sorts the original reference points by rho_St, but just the index, not the physical location
  */
void sort_all_refpoints_by_rho_index()
{
	int i,j,k;
	double temp, temp_index;
	int rho_sorted[factorial-adaptive_ref_points_inserted+factorial_inside];
	/*counter=0;
	for (i=0;i<factorial+factorial_inside;i++)
	{
		counter+=rho_St[i];
		printf("counter is %d\n",rho_St[i]);
	}*/
	for (i=0;i<factorial-adaptive_ref_points_inserted+factorial_inside;i++)
	{

		rho_sorted[i]=rho_St[i];
		printf("rho_St_sorted[%d] %d\n",i,rho_sorted[i]);
		sort_all_refpoint_index[i]=i;
	}
	for (i=0;i<factorial-adaptive_ref_points_inserted+factorial_inside;i++)
	{
		for (j=0;j<factorial-adaptive_ref_points_inserted+factorial_inside-i-1;j++)
		{
			if (rho_sorted[j]<rho_sorted[j+1])
			{
				temp=rho_sorted[j];
				temp_index=sort_all_refpoint_index[j];
				rho_sorted[j]=rho_sorted[j+1];
				sort_all_refpoint_index[j]=sort_all_refpoint_index[j+1];
				rho_sorted[j+1]=temp;
				sort_all_refpoint_index[j+1]=temp_index;
			}
		}
	}
	printf("\nSorting original Das and Dennis reference points by rho index\n________________________________________________\n");
	printf("rho_st is Sorted from highest to lowest, using sort_all_refpoint_index[i]\n");
	for (i=0;i<factorial-adaptive_ref_points_inserted+factorial_inside;i++)
	{
		printf("i %d,i_sort %d, rho_sort %d, rho_St_sort %d, rho_Fl_sort %d, rho_St %d, rho_Fl %d, rho %d\n",i,sort_all_refpoint_index[i],rho[sort_all_refpoint_index[i]],rho_St[sort_all_refpoint_index[i]],rho_Fl[sort_all_refpoint_index[i]],rho_St[i],rho_Fl[i],rho[i]);
	}
return;
}
/*This function finds the index of the useless refpoints, the usefull and useless numbers
  It considers whether rho_St nor rho_fl are equal to zero*/
void find_useless_usefull_refpoint_index()
{
	int i;
	useless_refpoint_number=0;
	usefull_refpoint_number=0;
	printf("adaptive_ref_points_inserted %d\n",adaptive_ref_points_inserted);
	for (i=0;i<factorial-adaptive_ref_points_inserted+factorial_inside;i++)
	{
		if (rho_St[i]==0)
		{
				useless_refpoint_index[useless_refpoint_number]=i;
				/*printf("useless_refpoint_index[%d] %d, rho_St %d\n",useless_refpoint_number,
				useless_refpoint_index[useless_refpoint_number],rho_St[i]);*/
				useless_refpoint_number++;

		}
		else if (rho_St[i]>0)
		{
				usefull_refpoint_index[usefull_refpoint_number]=i;
				/*printf("usefull_refpoint_index[%d] %d, rho_St %d\n",usefull_refpoint_number,
				usefull_refpoint_index[usefull_refpoint_number],rho_St[i]);*/
				usefull_refpoint_number++;
		}
	}
	printf("usefull_refpoint_number %d\n",usefull_refpoint_number);
	printf("useless_refpoint_number %d\n",useless_refpoint_number);
return;
}
void store_useless_refpoints (int adaptive_ref_point_number, int useless_ref_point_number)
{
	int i,k;
	printf("Storing useless refpoints temporally.\nThose vectors are restored to the next generation:\n");

	for (i=0;i<useless_ref_point_number;i++)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			temp_refpoints[k][i]=ref_points[k][useless_refpoint_index[i]];
			printf("%e\t",temp_refpoints[k][i]);
		}	
		temp_refpoints_pointer[i]=useless_refpoint_index[i];	
		printf("\t--->%dth\n",temp_refpoints_pointer[i]);
	}
return;
}
void load_useless_refpoints(int adaptive_ref_point_number, int useless_ref_point_number)
{
	int i,k;
	printf("Loading useless refpoints to the next generation\n");
	for (i=0;i<useless_ref_point_number;i++)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			ref_points[k][temp_refpoints_pointer[i]]=temp_refpoints[k][i];
			printf("%e\t",ref_points[k][i]);
		}
		printf("\n");
	}
return;
}
