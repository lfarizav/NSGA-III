/* Nond-domination based selection routines */
/* The Copyright belongs to Luis Felipe Ariza Vesga (lfarizav@unal.edu.co). You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for research purposes. All publications which use this code should acknowledge the author. Luis Felipe Ariza Vesga. 
A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems. March, 2019. */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>

# include "global.h"
# include "rand.h"

/* Routine to perform non-dominated sorting */
void fill_nondominated_sort (population *selection_pop, population *mixed_pop, population *new_pop, int generation)
{
    int flag;
    int i, j, k;
    int end;
    int front_size;
    int archieve_size;
    int rank=1;
    list *pool;
    list *elite;
    list *temp1, *temp2;
    pool = (list *)malloc(sizeof(list));
    elite = (list *)malloc(sizeof(list));
    front_size = 0;
    archieve_size=0;
    pool->index = -1;
    pool->parent = NULL;
    pool->child = NULL;
    elite->index = -1;
    elite->parent = NULL;
    elite->child = NULL;
    temp1 = pool;

    /*bubble_sorting_infeasible_population_index() function, stores infeasible solutions index
    and sorts them from higher to lower constrain violation*/
    bubble_sorting_infeasible_population_index(mixed_pop);

    /*printf("number_is_infeasible %d, number_is_feasible %d, 2*popsize %d\n",number_is_infeasible,number_is_feasible,2*popsize);
    printf("\nThere are 4 cases for fill nondominated sort\n");*/
    /*printf("******************************************************************************************************************\n");
    printf("Case 1-> All solutions are infeasible\n");
    printf("Case 2-> Number_is_feasible < popsize\n");
    printf("Case 3-> No constrains or Number_is_feasible=2*popsize\n");
    printf("Case 4-> Number_is_feasible >= popsize\n");
    printf("Case 1 and 2 do not associate reference points. Cases 3 and 4 do not consider infeasible solutions at all\n");
    printf("*****************************************************************************************************************\n\n");*/

    if (number_is_infeasible==2*popsize || (number_is_feasible<popsize && number_is_feasible>0))
    {
	if (number_is_infeasible==2*popsize)
	{
			/*printf("\nCase 1. All solutions are infeasible. Solutions sorted by minimum constrain violation values.\n\n");
			printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",number_is_feasible,number_is_infeasible,popsize);*/
			k=0;
			for (i=0;i<popsize;i++)
			{
				copy_ind (&mixed_pop->ind[infeasible_population_sorted_list_index[i]], &new_pop->ind[i]);
			}

			/*printf("Visualization all infeasible solutions of the new_pop sorted by constraint violations\n");*/
			/*for (k=0;k<popsize; k++)
		    	{
				display_pop_ind_obj(&(new_pop->ind[k]),k);
		   	}*/
			return;
	}
	else
	{
		if (number_is_feasible<popsize && number_is_feasible>0)
		{
			/*printf("\nCase 2: Number_is_feasible < popsize, number_is_feasible %d\n\n",number_is_feasible);
			printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",number_is_feasible,number_is_infeasible,popsize);*/
			for (i=0;i<number_is_feasible;i++)
			{
				copy_ind (&mixed_pop->ind[feasible_population_sorted_list_index[i]], &new_pop->ind[i]);
			}

			for (i=number_is_feasible,j=0;i<popsize;i++,j++)
			{
				copy_ind (&mixed_pop->ind[infeasible_population_sorted_list_index[j]], &new_pop->ind[i]);
			}

			/*printf("Visualization feasible solutions of the selection_pop number_is_feasible<=popsize\n");*/
			/*for (k=0;k<number_is_feasible; k++)
			{
				display_pop_ind_obj(&(new_pop->ind[k]),k);
			}*/
			/*printf("Visualization infeasible solutions of the selection_pop number_is_feasible<=popsize\n");*/
			/*for (k=number_is_feasible;k<popsize; k++)
			{
				display_pop_ind_obj(&(new_pop->ind[k]),k);
			}*/

			/*printf("Visualization of the new_pop infeasible/feasible solutions\n");*/
			/*for (k=0;k<popsize; k++)
			{
				display_pop_ind_obj(&(new_pop->ind[k]),k);
			}*/
			return;
		}
	}
    }
    else
    {
	
	if (number_is_feasible>=popsize)
	{
	    if (number_is_feasible==2*popsize)
	    {
		/*printf("\nCase 3: No constrains\n\n");
		printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",number_is_feasible,number_is_infeasible,popsize);*/
		for (i=0; i<2*popsize; i++)
		{
			insert (temp1,i);
			temp1 = temp1->child;
		}
	    }
	    else
	    {
		if (number_is_feasible>=popsize && number_is_feasible<2*popsize)
		{
			/*printf("\nCase 4: Number_is_feasible>=popsize\n\n");
			printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",
			number_is_feasible,number_is_infeasible,popsize);*/
			for (i=0; i<number_is_feasible; i++)
			{
				insert (temp1,feasible_population_sorted_list_index[i]);
				temp1 = temp1->child;
			}
		}
	    }
	    i=0;
	    do
	    {
		temp1 = pool->child;
		insert (elite, temp1->index);
		front_size = 1;
		temp2 = elite->child;
		temp1 = del (temp1);
		temp1 = temp1->child;
		do
		{
		    temp2 = elite->child;
		    if (temp1==NULL)
		    {
		        break;
		    }
		    do
		    {
		        end = 0;
		        flag = check_dominance (&(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]));
		        if (flag == 1)
		        {
			    insert (pool, temp2->index);
		            temp2 = del (temp2);
		            front_size--;
		            temp2 = temp2->child;
		        }
		        if (flag == 0)
		        {
		            temp2 = temp2->child;
		        }
		        if (flag == -1)
		        {
		            end = 1;
		        }
		    }
		    while (end!=1 && temp2!=NULL);
		    if (flag == 0 || flag == 1)
		    {
		        insert (elite, temp1->index);
		        front_size++;
		        temp1 = del (temp1);
		    }
		    temp1 = temp1->child;
		}
		while (temp1 != NULL);
		temp2 = elite->child;
		j=i;
		/*printf("archieve_size %d, front_size %d\n",archieve_size,front_size);*/
		if ( (archieve_size+front_size) <= popsize )
		{
		    do
		    {
		        copy_ind (&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
		        copy_ind (&mixed_pop->ind[temp2->index], &selection_pop->ind[i]);
		        new_pop->ind[i].rank = rank;
		        archieve_size+=1;
		        temp2 = temp2->child;
		        i+=1;
			if (archieve_size==popsize)
			{
				printf("Lucky generation, archieve_size==popsize, Nothing to do. Jump to the next generation\n");
				return;
			}
		    }
		    while (temp2 != NULL);
		    fronts[rank-1]=archieve_size;
		    /*printf("front %d has %d individuals\n",rank,fronts[rank-1]);
		    printf("archieve_size %d, front_size %d\n",archieve_size,front_size);*/
		    if (j==0)
		    {
			first_front=archieve_size;
			/*printf("first front is %d\n",first_front);*/
		    }
		    rank+=1;
		}
		else
		{
		    /*check 1*/
		    /*visualization of the selection_pop without the last front*/
		    /*printf("1. visualization of the selection_pop without the last front in the fillnds.c file\n");*/
		    /*for (k=0;k<archieve_size; k++)
	    	    {
			display_pop_ind_obj(&(selection_pop->ind[k]),k);
	   	    }*/
		    associated_reference_points_fill (selection_pop, mixed_pop, new_pop, front_size,archieve_size, elite, generation);
		    /*printf("archieve_size is %d, front_size is %d, i or count is %d, j is %d\n",archieve_size,front_size,i,j);*/

		    archieve_size = popsize;
		    for (j=i; j<popsize; j++)
		    {
		        new_pop->ind[j].rank = rank;
		    }
		}
		temp2 = elite->child;
		do
		{
		    temp2 = del (temp2);
		    temp2 = temp2->child;
		}
		while (elite->child !=NULL);
	    }
	    while (archieve_size < popsize);
	    while (pool!=NULL)
	    {
		temp1 = pool;
		pool = pool->child;
		free (temp1);
	    }
	    while (elite!=NULL)
	    {
		temp1 = elite;
		elite = elite->child;
		free (temp1);
	    }
	    return;
	}
    }
}

/* Routine to fill a population with individuals associated with reference points (Das and Dennis's) and (Deb and Jain)*/
void associated_reference_points_fill (population *selection_pop, population *mixed_pop, population *new_pop, int front_size, int archieve_size, list *elite, int generation)
{
    int j,k,l;
    int archieve_and_front_sizes;
    int pop_size;
    int member_number;
    int temp_rho_St_total;
    int temp_rho_Fl_total;

    /*initialization of selection variables*/
    archieve_and_front_sizes=variables_initialization(selection_pop, mixed_pop, new_pop, front_size, archieve_size, elite);/*ok*/
    if (archieve_size==popsize)
    {
	printf("Lucky generation, archieve_size==popsize\n");
	return;
    } 

    /*Load first front to pop_size*/
    pop_size=(archieve_size>0 && archieve_size<popsize)?archieve_size:archieve_and_front_sizes;
    /*printf("archieve_size is %d, front_size is %d\n",archieve_size,front_size);*/

    /* Routine to find min from functions*/
    for (j=0; j<pop_size; j++)
    {
        find_min_from_functions(&(selection_pop->ind[j]),j,1);/*ok*/
	find_max_from_functions(&(selection_pop->ind[j]),j,1);/*ok*/
    }

    /*visualization of min or ideal points*/
    /*for (j=0; j<nobj; j++)
    {
    	printf("4. scale_obj_min[%d][index=%d] %e\n",j,scale_obj_min_ref,scale_obj_min[j]);
    }*/
    /*visualization of max point*/
    /*for (j=0; j<nobj; j++)
    {
    	printf("4. scale_obj_max[%d][index=%d] %e\n",j,scale_obj_max_ref,scale_obj_max[j]);
    }*/
    
    /*Routine that substracts the zmin to solutions*/
    for (k=0; k<archieve_and_front_sizes; k++)
    {
        obj_minus_zmin(&(selection_pop->ind[k]));/*ok*/
    }

    /*Visualization of the selection_pop minus zmin in the fillnds.c file
    printf("5. visualization of the whole selection_pop minus zmin in the fillnds.c file\n");
    for (j=0;j<archieve_and_front_sizes; j++)
    {
	display_pop_ind_obj_minus_zmin(&(selection_pop->ind[j]),j);
    }*/

    find_extreme_points(selection_pop, pop_size);/*only consider the individuals in the first front*/
    /*normalize solutions*/
    for (l=0; l<archieve_and_front_sizes; l++)
    {
        normalized_objective_function(&(selection_pop->ind[l]));/*ok*/
    }

    /*7. visualization of the selection_pop normalized in the fillnds.c file
    printf("7. visualization of the whole selection_pop normalized in the fillnds.c file\n");
    for (l=0; l<archieve_and_front_sizes; l++)
    {
 	display_pop_ind_obj_normalized (&(selection_pop->ind[l]),l);
    }*/

    /*Associate normalized solutions with reference points*/
    for (l=0; l<archieve_and_front_sizes; l++)
    {
        associate(&selection_pop->ind[l],l, archieve_size,0,factorial+factorial_inside+last_gen_adaptive_refpoints_number);
    }
    /*printf("Printing rho after asociation refpoints ((factorial %d-adaptive_ref_points_inserted %d=%d)+factorial_inside %d +last_gen_adaptive_refpoints_number %d)\n",factorial,adaptive_ref_points_inserted,factorial-adaptive_ref_points_inserted,factorial_inside,last_gen_adaptive_refpoints_number);*/
    temp_rho_St_total=0;
    temp_rho_Fl_total=0;
    for (l=0; l<factorial+factorial_inside+last_gen_adaptive_refpoints_number; l++)
    {
	temp_rho_St_total+=rho_St[l];
	temp_rho_Fl_total+=rho_Fl[l];
	/*printf("rho_St[%d] %d, rho_Fl[%d] %d\n",l,rho_St[l],l,rho_Fl[l]);*/
    }
    /*printf("After association: temp_rho_St_total %d, temp_rho_Fl_total %d\n",temp_rho_St_total,temp_rho_Fl_total);*/
    
    /*Visualization of associated and distance reference
    for (l=0; l<archieve_and_front_sizes; l++)
    {
 	printf("%d\tassociatedref %d, distancetoassociatedref %e\n",l,
	selection_pop->ind[l].associatedref,selection_pop->ind[l].distancetoassociatedref);
    }*/

    /*Visualization of associated and distance reference
    for (l=0; l<archieve_and_front_sizes; l++)
    {
 	printf("%d\tassociatedref %d, distancetoassociatedref %e\n",l,
	selection_pop->ind[l].associatedref,selection_pop->ind[l].distancetoassociatedref);
    }*/

    /*Select solutions from last_front to complete new_pop*/
    member_number=niching(selection_pop, new_pop, front_size,archieve_size,0,factorial+factorial_inside+last_gen_adaptive_refpoints_number);
    /*printf("After niching\n");*/
    temp_rho_St_total=0;
    temp_rho_Fl_total=0;
    for (l=0; l<factorial+factorial_inside+last_gen_adaptive_refpoints_number; l++)
    {
	temp_rho_St_total+=rho_St[l];
	temp_rho_Fl_total+=rho_Fl[l];
	/*printf("rho_St[%d] %d, rho_Fl[%d] %d\t",l,rho_St[l],l,rho_Fl[l]);*/
	/*for (k=0;k<nobj;k++)
	{
		printf("%e\t",ref_points[k][l]);
	}
	printf("\n");*/
    }
    /*printf("After niching: temp_rho_St_total %d, temp_rho_Fl_total %d\n",temp_rho_St_total,temp_rho_Fl_total);
    printf("member_number %d, archieve_size %d,total %d, popsize %d\n",member_number,archieve_size,member_number+archieve_size,popsize);*/


    /*If adaptive nsga-iii is enabled, then,*/
    if (adaptive_nsga==1 || adaptive_nsga==2)
    {
	/*Add adaptive refpoints to improve Pareto Front distribution*/
	add_adaptive_refpoints_to_ref_points();
	last_gen_adaptive_refpoints_number=delete_adaptive_refpoints(archieve_size,front_size,selection_pop,new_pop,generation);
	/*display_refpoints ();*/
	/*initialization of rho's for next association*/
        /*for (j=0;j<factorial+factorial_inside;j++)
        {
		rho_St[j]=0;
		rho_Fl[j]=0;
		rho[j]=0;
        }*/

	/*Associate normalized solutions with reference and adaptive points*/
	/*for (l=0; l<archieve_and_front_sizes; l++)
	{
		associate(&(selection_pop->ind[l]),&(new_pop->ind[l]),l, archieve_size);
	}*/
	/*Preserve the rho, before niching change its values*/
        /*for (j=0;j<factorial+factorial_inside;j++)
        {
		if (archieve_size==0)
			rho_St_adaptive[j]=rho_Fl[j];
		else
			rho_St_adaptive[j]=rho_St[j];
        }*/

    }
    /*If adaptive nsga-iii is enabled, then,*/
    /*if (adaptive_nsga==1 || adaptive_nsga==2)
    {*/
    	/*Delete useless adaptive refpoints*/
    	/*last_gen_adaptive_refpoints_number=delete_adaptive_refpoints();

    }*/
}
int bubble_sorting_infeasible_population_index(population *poputation_sorted)
{
	int i,j,k;
	double temp, temp_index;
	double CV[2*popsize];

    	number_is_feasible = 0;
    	number_is_infeasible = 0;


	for (i=0;i<2*popsize;i++)
	{
		/*Checking feasibility of mixed solutions*/
		if (poputation_sorted->ind[i].constr_violation+poputation_sorted->ind[i].equality_constr_violation<0)
		{
			infeasible_population_sorted_list_index[number_is_infeasible]=i;
			CV[number_is_infeasible]=poputation_sorted->ind[i].constr_violation+poputation_sorted->ind[i].equality_constr_violation;
			poputation_sorted->ind[i].is_feasible=0;
			number_is_infeasible++;
			/*printf("%d,is infeasible %d\n",i,poputation_sorted->ind[i].is_feasible);*/
		}
		else
		{
			feasible_population_sorted_list_index[number_is_feasible]=i;
			poputation_sorted->ind[i].is_feasible=1;
			number_is_feasible++;
			/*printf("%d,is feasible %d, number_is_feasible %d\n",i,poputation_sorted->ind[i].is_feasible,number_is_feasible);*/
		}
	}
	if (number_is_feasible==2*popsize)
		return number_is_feasible;
	for (i=0;i<number_is_infeasible;i++)
	{
		for (j=0;j<number_is_infeasible-i-1;j++)
		{
			if (CV[j]<CV[j+1])
			{
				temp=CV[j];
				temp_index=infeasible_population_sorted_list_index[j];
				CV[j]=CV[j+1];
				infeasible_population_sorted_list_index[j]=infeasible_population_sorted_list_index[j+1];
				CV[j+1]=temp;
				infeasible_population_sorted_list_index[j+1]=temp_index;
			}
		}
	}
	/*printf("sorting infeasible population by violation constrains\n");*/
	/*for (i=0;i<number_is_infeasible;i++)
	{
		printf("%d,%e\n",i,poputation_sorted->ind[infeasible_population_sorted_list_index[i]].constr_violation+poputation_sorted->ind[infeasible_population_sorted_list_index[i]].equality_constr_violation);
	}*/
	/*printf("sorting feasible population by violation constrains\n");*/
	/*for (i=0;i<number_is_feasible;i++)
	{
		printf("%d,%e\n",i,poputation_sorted->ind[feasible_population_sorted_list_index[i]].constr_violation+poputation_sorted->ind[feasible_population_sorted_list_index[i]].equality_constr_violation);
	}*/
return number_is_feasible;
}
void feasible_population_index(population *poputation_sorted)
{
	int i,j,k;
	double temp, temp_index;
	double CV[2*popsize];
	k=0;
	for (i=0;i<2*popsize;i++)
	{
		if (poputation_sorted->ind[i].is_feasible)
		{
			feasible_population_sorted_list_index[k]=i;
			/*printf("feasible_population_index[i] %d\n",i);*/
			k++;
		}
	}
return;
}
void infeasible_population_index(population *poputation_sorted)
{
	int i,j,k;
	double temp, temp_index;
	double CV[2*popsize];
	k=0;
	for (i=0;i<2*popsize;i++)
	{
		if (!(poputation_sorted->ind[i].is_feasible))
		{
			infeasible_population_sorted_list_index[k]=i;
			/*printf("infeasible_population_index[i] %d\n",i);*/
			k++;
		}
	}
return;
}
