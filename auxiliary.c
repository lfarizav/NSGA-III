/* Some utility functions (not part of the algorithm) */
/* The Copyright belongs to Luis Felipe Ariza Vesga (lfarizav@unal.edu.co). You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for research purposes. All publications which use this code should acknowledge the author. Luis Felipe Ariza Vesga. 
A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Many-Objective Problems. March, 2019. */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to return the maximum of two variables */
double maximum (double a, double b)
{
    if (a>b)
    {
        return(a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
double minimum (double a, double b)
{
    if (a<b)
    {
        return (a);
    }
    return (b);
}
