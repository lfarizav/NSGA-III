/* Routine for evaluating population members  */
/* The Copyright belongs to Luis Felipe Ariza Vesga (lfarizav@unal.edu.co). You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for research purposes. All publications which use this code should acknowledge the author. Luis Felipe Ariza Vesga. 
A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems. March, 2019. */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine to evaluate objective function values and constraints for a population */
void evaluate_pop (population *pop)
{
    int i;
    /*printf("The following is the true pareto front of dtlz5\n");*/
    for (i=0; i<popsize; i++)
    {
        evaluate_ind (&(pop->ind[i]),i);
    }
    return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void evaluate_ind (individual *ind, int psize)
{
    int i,j,k;
    int normalized=0;
    test_problem (ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr, ind->equality_constr, normalized);
    /*for (i=0;i<nobj;i++)
    	printf("%e\t",ind->obj[i]);
    printf("\n");*/
    if (ncon==0 && neqcon==0)
    {
        ind->constr_violation = 0.0;
	ind->equality_constr_violation = 0.0;
    }
    else
    {
        ind->constr_violation = 0.0;
	ind->equality_constr_violation = 0.0;
        for (j=0; j<ncon; j++)
        {
            if (ind->constr[j]<0.0)
            {
                ind->constr_violation += ind->constr[j];
            }
        }
	if (neqcon>0)
	{
		for (k=0; k<neqcon; k++)
		{
		        ind->equality_constr_violation += abs(ind->equality_constr[k]);
		}
	}
	else
		ind->equality_constr_violation = 0;

    }
    return;
}

