/* Data initializtion routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include <string.h>
# include "global.h"

int variables_initialization(population *selection_pop, population *mixed_pop, population *new_pop, int front_size, int archieve_size, list *elite)
{
    int i,j;
    int archieve_and_front_sizes;
    int dist[2*popsize];
    list *temp;
    temp = elite->child;

    for (i=0;i<nobj;i++)
    {
	scale_obj_min[i]=DBL_MAX;
    }
    memset(scale_obj_max,0,nobj*sizeof(double));
    memset(scale_obj_min_ref,0,nobj*sizeof(int));
    memset(scale_obj_max_ref,0,nobj*sizeof(int));

    for (i=0;i<nobj;i++)
    {
        smin[i]=DBL_MAX;
        index_s[i]=0;
        a[i]=0;
        for (j=0;j<nobj;j++)
        {
            zmax[i][j]=0;
        }      
    }

    for (i=0;i<2*popsize;i++)
    {
	associatedfromlastfront_St[i]=0;
	associatedfromlastfront_Fl[i]=0;
	membertoadd[i]=0;
    }
    /*initialization of rho's*/
    for (j=0;j<factorial+factorial_inside;j++)
    {
	rho_St[j]=0;
	rho_Fl[j]=0;
	rho[j]=0;
    }
    for (j=0;j<factorial+factorial_inside;j++)
    {
	rho_St_adaptive[j]=0;
	rho_Fl_adaptive[j]=0;
	rho_adaptive[j]=0;
    }
    archieve_and_front_sizes=front_size+archieve_size;

    if (number_is_infeasible>2*popsize || number_is_feasible>2*popsize)
    {
	printf("error in eval_individual function, number_is_infeasible is %d number_is_feasible is %d, greater than %d\n",number_is_infeasible,number_is_feasible,2*popsize);
	exit(-1);
    }

    if (number_is_feasible>popsize && number_is_feasible<=2*popsize)
    {
	if (number_is_feasible==2*popsize)
		printf("Case 3: No constrains\n");
	else
		printf("Case 4: Number_is_feasible>popsize\n");
	for (j=archieve_size; j<archieve_and_front_sizes; j++)
	{
		dist[j] = temp->index;
		temp = temp->child;
		copy_ind(&mixed_pop->ind[dist[j]], &selection_pop->ind[j]);	
	}
	/*Visualization of selection pop*/
	printf("Visualization of the selection pop (archieved + front sizes)\n");
	for (j=0; j<archieve_and_front_sizes; j++)
	{
		display_pop_ind_obj (&(selection_pop->ind[j]),j);
	}
    }  
return archieve_and_front_sizes;  
}

void find_min_from_functions(individual *ind, int k, int population_type)
{
    int i;
    int normal_population=0;
    int minus_zmin_population=0;
    int normalized_population=0;

    if (population_type==1)
	normal_population=1;
    if (population_type==2)
	minus_zmin_population=1;
    if (population_type==3)
	normalized_population=1;
  
    for (i=0;i<nobj;i++)
    {
	if (normal_population)
	{
		if (ind->obj[i]<scale_obj_min[i])
		{
		    scale_obj_min[i]=ind->obj[i];
		    scale_obj_min_ref[i]=k;
		    /*printf("scale_obj_min[%d] = %e\n",i,scale_obj_min[i]);*/
		}
	}
	if (minus_zmin_population)
	{
		if (ind->obj_minus_zmin[i]<scale_obj_min[i])
		{
		    scale_obj_min[i]=ind->obj_minus_zmin[i];
		    scale_obj_min_ref[i]=k;
		    /*printf("scale_obj_min[%d] = %e\n",i,scale_obj_min[i]);*/
		}
	}
	if (normalized_population)
	{
		if (ind->obj_normalized[i]<scale_obj_min[i])
		{
		    scale_obj_min[i]=ind->obj_normalized[i];
		    scale_obj_min_ref[i]=k;
		    /*printf("scale_obj_min[%d] = %e\n",i,scale_obj_min[i]);*/
		}
	}

    }
    return;
}
void find_max_from_functions(individual *ind, int k, int population_type)
{
     int i;
    int normal_population=0;
    int minus_zmin_population=0;
    int normalized_population=0;

    if (population_type==1)
	normal_population=1;
    if (population_type==2)
	minus_zmin_population=1;
    if (population_type==3)
	normalized_population=1;
  
    for (i=0;i<nobj;i++)
    {
	if (normal_population)
	{
		if (ind->obj[i]>scale_obj_max[i])
		{
		    scale_obj_max[i]=ind->obj[i];
		    scale_obj_max_ref[i]=k;
		    /*printf("scale_obj_max[%d] = %e\n",i,scale_obj_max[i]);*/
		}
	}
	if (minus_zmin_population)
	{
		if (ind->obj_minus_zmin[i]>scale_obj_max[i])
		{
		    scale_obj_max[i]=ind->obj_minus_zmin[i];
		    scale_obj_max_ref[i]=k;
		    /*printf("scale_obj_max[%d] = %e\n",i,scale_obj_max[i]);*/
		}
	}
	if (normalized_population)
	{
		if (ind->obj_normalized[i]>scale_obj_max[i])
		{
		    scale_obj_max[i]=ind->obj_normalized[i];
		    scale_obj_max_ref[i]=k;
		    /*printf("scale_obj_max[%d] = %e\n",i,scale_obj_max[i]);*/
		}
	}

    }
    return;
}

void obj_minus_zmin (individual *ind)
{
    int i;
    for(i=0;i<nobj;i++)
    {
	ind->obj_minus_zmin[i]=ind->obj[i]-scale_obj_min[i];
    }
    return;
}
void associate(individual *normalizedind,individual *new_ind, int l, int archieve_size)
{
    int i,j,k;
    int index_dmin;
    double pd;
    double dmin;
    dmin=DBL_MAX;

    for (i=0;i<factorial+factorial_inside;i++)
    {
	pd=perpendicular_distance(normalizedind, i);
	if ( pd< dmin)
	{
		dmin = pd;
		index_dmin= i;
	}
    }
    normalizedind->associatedref=index_dmin;
    normalizedind->distancetoassociatedref=dmin;
    normalizedind->w=num_div_den[index_dmin];
    rho[index_dmin]+=1;
    if (l<archieve_size)
    {
	rho_St[index_dmin]+=1;	
	
	/*The follow information is load from selection_pop to new_pop. Just for visualization purposes*/
	new_ind->associatedref=index_dmin;
    	new_ind->distancetoassociatedref=dmin;
    	new_ind->w=num_div_den[index_dmin];
	for (k=0;k<nobj;k++)
	{
		new_ind->obj_normalized[k]=normalizedind->obj_normalized[k];
		new_ind->obj_minus_zmin[k]=normalizedind->obj_minus_zmin[k];
	}
	/**********************************************************************************************/
    }
    else
	rho_Fl[index_dmin]+=1;
return;
}

double perpendicular_distance(individual *normalizedind, int l)
{ 
	double numerator=0.0, denominator=0.0, d;
	int k;

	for (k=0;k<nobj;k++)
	{
		numerator += ref_points[k][l]*normalizedind->obj_normalized[k];
		denominator += pow(ref_points[k][l],2);
	}
	/*num_div_den is used to visualize the perpendicular vector from one solution to the reference line*/
	num_div_den[l]=numerator/denominator;
	d=0;
	for (k=0;k<nobj;k++)
	{
		d+=pow(num_div_den[l]*ref_points[k][l]-normalizedind->obj_normalized[k],2);
	}
return sqrt(d);
}

void niching (population *selection_pop,population *new_pop, int front_size, int archieve_size)
{
    int min_rho_index, min_rho_St_index, min_rho_Fl_index;
    int numberofindexassociatedtorefpoint_j;
    int i,l;
    int rand;
    double min_ddj;
    int associatedfromlastfront_index;
    int new_member_index_from_Fl;
    int membernumber;
    int index_St_and_Fl=0;
    int index_Fl=0;
    membernumber=0;
    do{
	    min_rho_index=0;
	    min_rho_St=DBL_MAX;
	    min_rho_Fl=DBL_MAX;
	    min_rho=DBL_MAX;
	    /*Finds minimun rho*/
	    for (i=0; i<factorial+factorial_inside; i++)
	    {
	        if (rho_Fl[i]<min_rho_Fl)
		{
		    min_rho_Fl=rho_Fl[i];
		}
	        if (rho[i]<min_rho)
		{
		    min_rho=rho[i];
		    min_rho_index=i;
		}
		/*If rho or rho_Fl are zero, then discard the reference point associated*/
		if (rho[i]==0 || rho_Fl[i]==0)/*Second scenario: min_rho_St=0, min_rho_Fl=0. Excluded from further consideration for the current generation*/
		{
		    rho_St[i]=DBL_MAX;
		    rho_Fl[i]=DBL_MAX;
		    rho[i]   =DBL_MAX;
		}
	    }

	    if (min_rho==0 || min_rho_Fl==0)
	    {
		/*Visualization without excluding reference points with rho==0
	    	for (l=0; l<factorial+factorial_inside; l++)
	    	{
	 		printf("rho[%d] is %d, rho_Fl[%d] is %d\n",l,rho[l],l,rho_Fl[l]);

	    	}
	    	printf("min_rho_Fl is %d,min_rho is %d, min_rho_index is %d \n",min_rho_Fl,min_rho,min_rho_index);*/
		continue;
  	    }
	    /*Visualizaton excluding reference points with rho==0
	    for (l=0; l<factorial+factorial_inside; l++)
	    {
		if(rho[l]<100000000)
		 	printf("rho_St[%d] %d, rho_Fl[%d] %d, rho[%d] %d\n",l,rho_St[l],l,rho_Fl[l],l,rho[l]);
	    }*/
	    /*Finds the index for all reference points with the same min_rho*/
	    index_St_and_Fl=0;
	    for (i=0; i<factorial+factorial_inside; i++)
	    {
		if (rho[i]==min_rho)
		{
			ref_points_min_rho[index_St_and_Fl]=i;
			index_St_and_Fl++;
		}
	    }
	    /*If index is >1, then select randomly a index which has the same min_rho*/
	    if (index_St_and_Fl>1)
	    {
		/*printf("rand is %d, ref_points_min_rho[%d] is %d\n",rand,rand,ref_points_min_rho[rand]);*/
		rand=rnd(0,index_St_and_Fl-1);
	    	min_rho_index=ref_points_min_rho[rand];
	    }
		
	    numberofindexassociatedtorefpoint_j=0;
	    for (l=archieve_size; l<archieve_size+front_size; l++)
	    {
	    	if (associated_from_last_front(&(selection_pop->ind[l]),l,min_rho_index,numberofindexassociatedtorefpoint_j)==1)
		{
			numberofindexassociatedtorefpoint_j++;
		}
	    }

	    /*printf("numberofindexassociatedtorefpoint_j is %d\n",numberofindexassociatedtorefpoint_j);*/
	    new_member_index_from_Fl=0;    
	    min_ddj=DBL_MAX;
	    for (l=0; l<numberofindexassociatedtorefpoint_j; l++)
	    {
			/*printf("d[%d][%d] is %e, distancetoassociatedref is %e, min_ddj is %e\n",associatedfromlastfront_Fl[l],min_rho_index,d[associatedfromlastfront_Fl[l]][min_rho_index],selection_pop->ind[associatedfromlastfront_Fl[l]].distancetoassociatedref,min_ddj);*/
			if (selection_pop->ind[associatedfromlastfront_Fl[l]].distancetoassociatedref<min_ddj)
			{
					min_ddj=selection_pop->ind[associatedfromlastfront_Fl[l]].distancetoassociatedref;
					new_member_index_from_Fl=associatedfromlastfront_Fl[l];/*index in the selection pop*/
			}	
	    }
	    /*printf("The minimun distance is %e, the reference is %d\n",min_ddj,new_member_index_from_Fl);*/
	    selection_pop->ind[new_member_index_from_Fl].associatedref=-1;
	    membertoadd[membernumber]=new_member_index_from_Fl;
	    membernumber++;

	    if (numberofindexassociatedtorefpoint_j==1)	
	    {
		    rho_St[min_rho_index]=DBL_MAX;
		    rho_Fl[min_rho_index]=DBL_MAX;
		    rho[min_rho_index]=DBL_MAX;
	    }
	    else if (numberofindexassociatedtorefpoint_j>1)	
	    {
	    rho_Fl[min_rho_index]+=1;
	    rho[min_rho_index]   +=1;
	    }

	    /*printf("associatedref is %d, min_rho_index is %d\n",selection_pop->ind[new_member_index_from_Fl].associatedref,min_rho_index);*/
	    /*selection_pop->ind[new_member_index_from_Fl].associatedref=DBL_MAX;
	    printf("associatedref is %d\n",selection_pop->ind[new_member_index_from_Fl].associatedref);*/

	    if (archieve_size+membernumber>=popsize)
		break;
    }while (1);

    /*printf("member number is %d\n",membernumber);
    for (l=0; l<popsize; l++)
    {
	display_pop_ind_obj(&(new_pop->ind[l]),l);
    }*/
    for (l=0; l<membernumber; l++)
    {
	copy_ind(&selection_pop->ind[membertoadd[l]], &new_pop->ind[(archieve_size==0)?l:archieve_size+l]);
    }
}
int associated_from_last_front(individual *normalizedind, int l, int index, int associatedfromlastfront_index)
{
	if (normalizedind->associatedref==index)
	{
		associatedfromlastfront_Fl[associatedfromlastfront_index]=l;/*ind associated from last front*/
		return 1;	
	}
	else{
		return 0;
	}
}
void normalized_objective_function (individual *ind)/*ok*/
{
    int i;
    for (i=0;i<nobj;i++)
    {
	if (a[i]>0.0000000001)
		ind->obj_normalized[i]=ind->obj_minus_zmin[i]/a[i];
	else
		ind->obj_normalized[i]=ind->obj_minus_zmin[i]/0.0000000001;
    }
    return;

}
void normalized_objective_function_simple (individual *ind)
{
    int i;
    for (i=0;i<nobj;i++)
    {
	ind->obj_normalized[i]=ind->obj_minus_zmin[i]/(scale_obj_max[i]-scale_obj_min[i]);
    }
    return;
}

void find_a ()/*ok*/
{
    
    int i,j;
    double zmax_matrix[25][25];
    double d=0.0;
    /*check 6*/
    /*6. matrix zmax is*/
    printf("6. zmax matrix is: \n");
    for (i=0;i<nobj;i++)
    {
        for(j=0;j<nobj;j++)
        {
            zmax_matrix[i][j]=zmax[i][j];
            printf("\t%e\t",zmax_matrix[i][j]);
        }
        printf("\n");
    }
    d = determinant(zmax_matrix, nobj);
    /*printf("determinant is %e\n",d);*/
    if (d == 0)
        printf("\nInverse matrix impossible. Singular Matrix\n");
    else
        cofactor(zmax_matrix, nobj);
    return;
}

int is_zmax_duplicated ()/*ok*/
{
	/* Check whether there are duplicate extreme points.
	This might happen but the original paper does not mention how to deal with it.*/
	int i,j,k;
	int is_duplicated = 0;
	for (i=0; !is_duplicated && i<nobj; i++)
	{
		for (j=i+1; !is_duplicated && j<nobj; j++)
		{
			int count_duplicate_per_dim=0;
			for (k=0; k<nobj; k++)
			{
				if (zmax[i][k]==zmax[j][k])
					count_duplicate_per_dim+=1;
			}
        		/*printf("count_duplicate_per_dim is %d\n",count_duplicate_per_dim);*/
			if (count_duplicate_per_dim==nobj)
				is_duplicated = 1;
		}
	}

	return is_duplicated;
}
void construct_hyperplane (population *selection_pop, int pop_size)/*ok*/
{
	int i,j,k;
	int duplicated, negative_intercept;
	negative_intercept=0;
	duplicated=is_zmax_duplicated ();
	if (!duplicated)
		printf("duplicate %d\n",duplicated);
	{
		find_a();
		for (i=0;i<nobj;i++)
    		{
			if (a[i]<0)
			{
				negative_intercept = 1;
				break;
			}
		}
	}
	if (duplicated || negative_intercept)
	{
		printf("duplicate %d or negative_intercept %d\n",duplicated,negative_intercept);
	    	/*for (i=0;i<nobj;i++)
	    	{
			scale_obj_max[i]=0;
	    	}

		for (j=0; j<pop_size; j++)
    		{
        		find_max_from_functions(&(selection_pop->ind[j]),j,1);
    		}*/

		printf("When duplicated or negative intercept, the a vector is:\n");
		for (i=0;i<nobj;i++)
	    	{
			a[i]=selection_pop->ind[scale_obj_max_ref[i]].obj[i];
			printf("a[%d]= %e\n",i,a[i]);
	    	}
		printf("\n");
	}
}
double achievement_scalarization_function (individual *ind_minus_zmin,int i)/*ok*/
{
    int k;
    double temp=0;
    for (k=0;k<nobj;k++)
    {
	if (ind_minus_zmin->obj_minus_zmin[k]/w_scalarizing_vector[k]>temp)
	{
	   temp=ind_minus_zmin->obj_minus_zmin[k]/w_scalarizing_vector[k];
	}	
    }
    printf("%d\t",i);
    for (k=0;k<nobj;k++)
    {
        printf("%e\t",ind_minus_zmin->obj_minus_zmin[k]/w_scalarizing_vector[k]);
    }
    printf("max --> %e\n", temp);
    return temp;
}
void find_extreme_points(population *selection_pop_minus_zmin, int pop_size)/*ok*/
{
    int i,j,k,l,m;
    double temp_s_min;
    int temp_s_min_index;
    double s;

	for (i=0;i<nobj;i++)
	{
		get_scalarizing_vector(i);
		temp_s_min=DBL_MAX;
		for (j=0;j<pop_size;j++)
		{
		    s=achievement_scalarization_function(&(selection_pop_minus_zmin->ind[j]), j);
		    if (s<temp_s_min)
		    {
		        temp_s_min=s;
		        temp_s_min_index=j;
		    }
		}
		for (m=0;m<nobj;m++)
		   {
		       zmax[i][m]=selection_pop_minus_zmin->ind[temp_s_min_index].obj_minus_zmin[m];
		       printf("%e\t",selection_pop_minus_zmin->ind[temp_s_min_index].obj_minus_zmin[m]);
		   }
		    printf("\n");
	}
	construct_hyperplane (selection_pop_minus_zmin,pop_size);
    return;
}
void get_scalarizing_vector(int j)/*ok*/
{
    double epsilon = 0.0000000001;
    int i;
    /*printf("The scalarization vector is:\n");*/
    for (i=0;i<nobj;i++)
    {
	if (i==j)
		w_scalarizing_vector[i]=1;
	else
		w_scalarizing_vector[i]=epsilon;
	/*printf("w_scalarizing_vector[%d] %e\n",i,w_scalarizing_vector[i]);*/
    }
    return;
}

/*For calculating Determinant of the Matrix */
double determinant(double zmax_matrix[25][25], int nobj)/*ok*/
{
    double s = 1, det = 0, b[25][25];
    int i, j, m, n, c;
    if (nobj == 1)
    {
        return (zmax_matrix[0][0]);
    }
    else
    {
        det = 0;
        for (c = 0; c < nobj; c++)
        {
            m = 0;
            n = 0;
            for (i = 0;i < nobj; i++)
            {
                for (j = 0 ;j < nobj; j++)
                {
                    b[i][j] = 0;
                    if (i != 0 && j != c)
                    {
                        b[m][n] = zmax_matrix[i][j];
                        if (n < (nobj - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (zmax_matrix[0][c] * determinant(b, nobj - 1));
            s = -1 * s;
        }
    }
    
    return (det);
}
void cofactor(double num[25][25], int f)/*ok*/
{
    double b[25][25], fac[25][25];
    int p, q, m, n, i, j;
    for (q = 0;q < f; q++)
    {
        for (p = 0;p < f; p++)
        {
            m = 0;
            n = 0;
            for (i = 0;i < f; i++)
            {
                for (j = 0;j < f; j++)
                {
                    if (i != q && j != p)
                    {
                        b[m][n] = num[i][j];
                        if (n < (f - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
        }
    }
    transpose(num, fac, f);
}
/*Finding transpose of matrix*/
void transpose(double num[25][25], double fac[25][25], int r)/*ok*/
{
    int i, j;
    float b[25][25], inverse[25][25], d;
    
    for (i = 0;i < r; i++)
    {
        for (j = 0;j < r; j++)
        {
            b[i][j] = fac[j][i];
        }
    }
    d = determinant(num, r);
    for (i = 0;i < r; i++)
    {
        for (j = 0;j < r; j++)
        {
            inverse[i][j] = b[i][j] / d;
        }
    }
    printf("6.  The inverse zmax matrix is : \n");
    
    for (i = 0;i < r; i++)
    {
        for (j = 0;j < r; j++)
        {
            printf("\t%f", inverse[i][j]);
            a[i]+=inverse[j][i];
        }
        printf("\n");
    }
    for (j = 0;j < r; j++)
    {
        a[j]=1/a[j];
        printf("a[%d] is %e\n",j,a[j]);
    }
}

