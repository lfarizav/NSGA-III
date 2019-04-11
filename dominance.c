/* Domination checking routines */
/* The Copyright belongs to Luis Felipe Ariza Vesga (lfarizav@unal.edu.co). You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for research purposes. All publications which use this code should acknowledge the author. Luis Felipe Ariza Vesga. 
A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Many-Objective Problems. March, 2019. */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine for usual non-domination checking
   It will return the following values
   1 if a dominates b
   -1 if b dominates a
   0 if both a and b are non-dominated */

int check_dominance (individual *a, individual *b)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;

	if (a->is_feasible>1 && b->is_feasible>1 && ncon>0 && neqcon>0)
	{
		printf("error. Please check evaluate_ind function. a->is_feasible>1 or b->is_feasible>1 or both\n");
		printf("a->is_feasible is %d, b->is_feasible is %d, ncon is %d, neqcon is %d\n",a->is_feasible,b->is_feasible,ncon,neqcon);
		exit(-1);
    	}
 	if (!(a->is_feasible>1 && b->is_feasible>1))
	{
	    if (a->is_feasible==1 && b->is_feasible==0)
	    {
			return (1);
	    }
	    else
	    {
		if (b->is_feasible==1 && a->is_feasible==0)
		{
			return (-1);
		}
	    }
	}

	if ((a->constr_violation+a->equality_constr_violation)<0 && (b->constr_violation+b->equality_constr_violation)<0)
	{

			if (a->constr_violation+a->equality_constr_violation > b->constr_violation+b->equality_constr_violation)
			{
			    return (1);
			    
			}
			else
			{
			    if (a->constr_violation+a->equality_constr_violation < b->constr_violation+b->equality_constr_violation)
			    {
				return (-1);
			    }
			    else
			    {
				return (0);
			    }
			}
	}
	else
	{
			if (a->constr_violation+a->equality_constr_violation < 0 && b->constr_violation+b->equality_constr_violation == 0)
			{
			    return (-1);
			    
			}
			else
			{
			    if (a->constr_violation+a->equality_constr_violation == 0 && b->constr_violation+b->equality_constr_violation <0)
			    {
				return (1);
			    }
			    else
			    {
				for (i=0; i<nobj; i++)
				{
				    if (a->obj[i] < b->obj[i])
				    {
				        flag1 = 1;

				    }
				    else
				    {
				        if (a->obj[i] > b->obj[i])
				        {
				            flag2 = 1;
				        }
				    }
				}
				if (flag1==1 && flag2==0)
				{
				    return (1);
				}
				else
				{
				    if (flag1==0 && flag2==1)
				    {
				        return (-1);
				    }
				    else
				    {
				        return (0);
				    }
				}
			    }
			}
	}
}
