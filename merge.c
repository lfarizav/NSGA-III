/* Routine for mergeing two populations */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine to merge two populations into one */
void merge(population *pop1, population *pop2, population *pop3)
{
    int i, k;

    for (i=0; i<popsize; i++)
    {
        copy_ind (&(pop1->ind[i]), &(pop3->ind[i]));
    }
    for (i=0, k=popsize; i<popsize; i++, k++)
    {
        copy_ind (&(pop2->ind[i]), &(pop3->ind[k]));
    }

    return;
}

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
void copy_ind (individual *ind1, individual *ind2)
{
    int i, j;
    ind2->rank = ind1->rank;
    ind2->constr_violation = ind1->constr_violation;
    ind2->equality_constr_violation = ind1->equality_constr_violation;
    ind2->associatedref = ind1->associatedref;
    ind2->distancetoassociatedref = ind1->distancetoassociatedref;
    ind2->is_feasible = ind1->is_feasible;
    ind2->w = ind1->w;

    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            ind2->xreal[i] = ind1->xreal[i];
        }
    }
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            ind2->xbin[i] = ind1->xbin[i];
            for (j=0; j<nbits[i]; j++)
            {
                ind2->gene[i][j] = ind1->gene[i][j];
            }
        }
    }
    for (i=0; i<nobj; i++)
    {
        ind2->obj[i] = ind1->obj[i];
	ind2->obj_minus_zmin[i] = ind1->obj_minus_zmin[i];
	ind2->obj_normalized[i] = ind1->obj_normalized[i];
    }
    if (ncon!=0)
    {
        for (i=0; i<ncon; i++)
        {
            ind2->constr[i] = ind1->constr[i];
        }
    }
    if (neqcon!=0)
    {
        for (i=0; i<neqcon; i++)
        {
            ind2->equality_constr[i] = ind1->equality_constr[i];
        }
    }

    return;
}
