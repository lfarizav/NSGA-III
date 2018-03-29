/* Data initializtion routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include <string.h>

# include "global.h"

/* Function to generate reference points*/
/*step = nobj,1/(double)(numberpointperdim-1),factorial*/
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
		return 0;
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
				ref_points[k][count+factorial]=0.5*step*index[k]+0.5/nobj;
				temp+=index[k];
			}
			ref_points[k][count+factorial]=0.5*(1-step*(temp))+0.5/nobj;
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
		return 0;
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
				temp+=index[k];
			}
			minimum_amount_ref_points[k][count]=(1/(double)(numberpointperdim-1))*(1-step*(temp));
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
	int i,j,k,l,m;
	int crowded_rho=1;
	int temp=0;
	int temp_j=last_gen_adaptive_refpoints_number;

	/*If last_gen_adaptive_refpoints_number is higher than factorial + factorial_inside, then do nothing */       
	if (temp_j>=factorial+factorial_inside)
		return temp_j;

	printf("Last generation, adaptive_refpoints starts from %d\n",temp_j);
	printf("Visualization of last generation adaptive refpoints\n");
	for (i=0;i<temp_j;i++)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			printf("%e\t",adaptive_refpoints[k][i]);
		}
		printf("\n");
	}
	int temp_nobj;
	int temp_nobj_adaptive;
	int temp_outside_first_quadrant;

	/*Organize refpoints from highest to lowest rho*/
    	sort_all_refpoints_by_rho_index();

	/*find_useless_refpoints_index finds index for useless_refpoints and initialize usefull_refpoint_number*/
        find_useless_usefull_refpoint_index();
	printf("The number of useless points is %d, the number of usefull points is %d\n", useless_refpoint_number, usefull_refpoint_number);

	/*Just check if rho>2 for usefull refpoints*/
	for (i=0;i<usefull_refpoint_number;i++)
	{
		printf("rho[%d] %d\n",sort_all_refpoint_index[i],rho[sort_all_refpoint_index[i]]);
		for (j=0;j<factorial_adaptive;j++)
		{
		    /*That is true because refpoints are sorted from highest to lowest values. So, 
	    	    it requires to do only usefull_refpoints_number times*/
		    if (rho[sort_all_refpoint_index[i]]>crowded_rho)
		    {
		      if (adaptive_nsga)
		      {
			for (k=0;k<nobj;k++)
			{
				/*Translation of minimum amount refpoints to the reference point which rho>1*/
				adaptive_refpoints[k][temp_j]=minimum_amount_ref_points[k][j]+
							ref_points[k][sort_all_refpoint_index[i]]-minimum_amount_ref_points[nobj-1][0]/nobj;
				/*printf("adap_refp[%d][%d] %e =\t",k,temp_j,adaptive_refpoints[k][temp_j]);
				printf("min_refp[%d][%d] %e +\t",k,temp_j,minimum_amount_ref_points[k][j]);
				printf("refp[%d][%d] %e -\t",k,temp_j,ref_points[k][sort_all_refpoint_index[i]]);
				printf("min_refp[nobj-1][0] %e\t",minimum_amount_ref_points[nobj-1][0]/nobj);
				printf("rho[%d] %d\n",i,rho[sort_all_refpoint_index[i]]);*/

			}
		      }
		      else if (adaptive_nsga==2)
		      {
			double translation[nobj];
			for (k=0;k<nobj;k++)
			{
				translation[k]=ref_points[k][sort_all_refpoint_index[i]]-minimum_amount_ref_points[k][j];
			}
			for (k=0;k<nobj;k++)
			{
				if (k!=j)
				{
					for (l=0;l<nobj;l++)
					{	
						adaptive_refpoints[l][temp_j]=minimum_amount_ref_points[l][k]+translation[l];
					}
				}
			}
		      }
			/*Verification:reference points adaptive neither outside the first quadrant nor repeated*/
			/*It either repeated of outside the first quadrant, assign the vector 0*/
			for (l=0;l<factorial+factorial_inside;l++)
			{
				temp_nobj=0;
				temp_outside_first_quadrant=0;
				for (k=0;k<nobj;k++)
				{
					if (adaptive_refpoints[k][temp_j]==ref_points[k][l])
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
					l=factorial+factorial_inside;
				}
			}
			/*Be sure there is not repeated adaptive refpoints, before add*/
			for (m=0;m<last_gen_adaptive_refpoints_number;m++)
			{
				temp_nobj_adaptive=0;
				for (k=0;k<nobj;k++)
				{
					if (adaptive_refpoints[k][temp_j]==adaptive_refpoints[k][m])
						temp_nobj_adaptive++;
				}
				if (temp_nobj_adaptive==nobj)
				{
					m=last_gen_adaptive_refpoints_number;
				}
			}	
			if (!(temp_nobj==nobj || temp_outside_first_quadrant || temp_nobj_adaptive==nobj))
			{
				printf("adaptive_refpoint [%d] added, after doing 3 checks\n",temp_j);
				temp_j++;
				/*The following if, restricts the size of adaptive_refpoint_number 
				to be less than factorial+factorial_inside*/
				if (temp_j>=factorial+factorial_inside)
				{
					i=usefull_refpoint_number;
				}
			}
		    }
		    else
		    {
			/*rho <2, then finish for and go to return*/
			i=usefull_refpoint_number;
		    }
		}
	}
	printf("The number of adaptive refpoints neither repeated nor outside the first quadrant is %d\n",temp_j-1);
	printf("Vacancies to fill, instead of original reference points with rho = 0 %d\n",useless_refpoint_number);
return (temp_j==0)?0:temp_j-1;
}

void add_adaptive_refpoints_to_ref_points()
{
	int i,j,k;
	int temp;
        temp_refpoints = (double **)malloc(nobj*sizeof(double *));
        for (i=0;i<nobj;i++)
        {
		temp_refpoints[i] = (double *)malloc(((nobj)*(factorial+factorial_inside))*sizeof(double));
        }
        /*create new adaptive reference points */
	/*Here adaptive_refpoint_number is returned*/
   	adaptive_refpoint_number=create_adaptive_refpoints();
        /*replace useless refpoints with adaptive refpoints*/

	for(i=0;i<adaptive_refpoint_number;i++)
	{
		printf("%d, rho[%d] %d\n",i,useless_refpoint_index[i],rho[useless_refpoint_index[i]]);
	}
	/*store refpoints replaced by adaptive refpoints*/
	store_useless_refpoints (adaptive_refpoint_number,useless_refpoint_number);
	for(i=0;i<adaptive_refpoint_number && i<useless_refpoint_number;i++)
	{
		for (k=0;k<nobj;k++)
		{
			ref_points[k][useless_refpoint_index[i]]=adaptive_refpoints[k][i];
			printf("adaptive_refpoint [%d][%d] %e, added when rho of refpoint is %d\n",k,i,adaptive_refpoints[k][i],rho[useless_refpoint_index[i]]);
		}
	}
return;
}
int delete_adaptive_refpoints()
{
	int i,j,k;
	int temp_adaptive_number=0;
	/*find useless and usefull refpoints
	find_useless_usefull_refpoint_index();*/
	for (i=0;i<adaptive_refpoint_number&& i<useless_refpoint_number;i++)
	{
		if (rho_adaptive[useless_refpoint_index[i]]>1)
		{
			for (k=0;k<nobj;k++)
			{
				adaptive_refpoints[k][temp_adaptive_number]=ref_points[k][useless_refpoint_index[i]];	
			}
			printf(" adaptive_refpoint [%d] added to next gen, rho_adaptive %d\n",temp_adaptive_number,rho_adaptive[useless_refpoint_index[i]]);
			temp_adaptive_number++;
		}
		else
		{
			printf("adaptive_refpoint [%d] deleted, rho_adaptive %d\n",i,rho_adaptive[useless_refpoint_index[i]]);
		}
	}
	printf("Visualization of last generation adaptive refpoints, before next generation\n");
	for (i=0;i<temp_adaptive_number;i++)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			printf("%e\t",adaptive_refpoints[k][i]);
		}
		printf("\n");
	}
return temp_adaptive_number;
}
/*This function sorts the reference points by rho, but just the index, not the physical location
  */
void sort_all_refpoints_by_rho_index()
{
	int i,j,k;
	double temp, temp_index;
	int rho_sorted[factorial+factorial_inside];

	for (i=0;i<factorial+factorial_inside;i++)
	{
		rho_sorted[i]=rho[i];
		sort_all_refpoint_index[i]=i;
	}
	for (i=0;i<factorial+factorial_inside;i++)
	{
		for (j=0;j<factorial+factorial_inside-i-1;j++)
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
	printf("rho is ordered from highest to lowest, using sort_all_refpoint_index[i]\n");
	for (i=0;i<factorial+factorial_inside;i++)
	{
		printf("i %d, rho %d\n",sort_all_refpoint_index[i],rho[sort_all_refpoint_index[i]]);
	}
return;
}
/*This function finds the index of the useless refpoints, the usefull and useless numbers
  It considers whether rho nor rho_fl are equal to zero*/
void find_useless_usefull_refpoint_index()
{
	int i;
	useless_refpoint_number=0;
	usefull_refpoint_number=0;
	for (i=0;i<factorial+factorial_inside;i++)
	{
		if (rho[i]==0)
		{
				useless_refpoint_index[useless_refpoint_number]=i;
				printf("useless_refpoint_index[%d] %d, rho %d\n",useless_refpoint_number,
				useless_refpoint_index[useless_refpoint_number],rho[i]);
				useless_refpoint_number++;

		}
		else
		{
			if (rho[i]>0)
			{
				usefull_refpoint_index[usefull_refpoint_number]=i;
				printf("usefull_refpoint_index[%d] %d, rho %d\n",usefull_refpoint_number,
				usefull_refpoint_index[usefull_refpoint_number],rho[i]);
				usefull_refpoint_number++;
			}
		}
	}
	printf("usefull_refpoint_number %d\n",usefull_refpoint_number);
	printf("useless_refpoint_number %d\n",useless_refpoint_number);
return;
}
void store_useless_refpoints (int adaptive_ref_point_number, int useless_ref_point_number)
{
	int i,k;
	printf("Storing useless refpoints temporally\n");
 	printf("test %e\n",ref_points[k][useless_refpoint_index[i]]);

	for (i=0;i<adaptive_ref_point_number && i<useless_ref_point_number;i++)
	{
		printf("%d\t",i);
		for (k=0;k<nobj;k++)
		{
			temp_refpoints[k][i]=ref_points[k][useless_refpoint_index[i]];
			printf("%e\t",temp_refpoints[k][i]);
		}		
		printf("\n");
	}
return;
}
void load_useless_refpoints(int adaptive_ref_point_number, int useless_ref_point_number)
{
	int i,k;
	printf("Loading useless refpoints to the next generation\n");
	for (i=0;i<adaptive_ref_point_number && i<useless_ref_point_number;i++)
	{
		for (k=0;k<nobj;k++)
		{
			ref_points[k][useless_refpoint_index[i]]=temp_refpoints[k][i];
			printf("%e\t",ref_points[k][i]);
		}
		printf("\n");
	}
return;
}
