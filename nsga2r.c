/* NSGA-II routine (implementation of the 'main' function) */
/* The Copyright belongs to Luis Felipe Ariza Vesga (lfarizav@unal.edu.co). You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for research purposes. All publications which use this code should acknowledge the author. Luis Felipe Ariza Vesga. 
A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems. March, 2019. */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>
# include <float.h>
# include <time.h>

# include "global.h"
# include "rand.h"

int nreal;
int nbin;
int nobj;
int ncon;
int neqcon;
int popsize;
int adaptive_nsga;
double *a_last_gen;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1;
int obj2;
int obj3;
int obj4;
int obj5;
int angle1;
int angle2;
double *scale_obj_min;
double *scale_obj_max;
int *scale_obj_min_ref;
int *scale_obj_max_ref;
int dtlz;
double *a;
double *smin;
int *index_s;
double **zmax;
long factorial;
long factorial_inside;
long factorial_adaptive;
int numberpointperdim;
int numberpointperdim_inside;
int *rho;
int *rho_St;
int *rho_Fl;
int *last_rho_St;
int *last_generation_associated_rho_St;
int min_rho_St;
int min_rho_Fl;
int min_rho;
int min_rho_Fl_count_index;
int min_rho_count_index;
int counter;
int *associatedfromlastfront_St;
int *associatedfromlastfront_Fl;
int *membertoadd;
int *ref_points_min_rho;
int *ref_points_min_rho_Fl;
double *w_scalarizing_vector;
int first_front;
int *fronts;
int *dist_lf;
double **ref_points;
double **adaptive_ref_points_settled;
int *adaptive_ref_points_settled_number;
int elegible_adaptive_ref_points_to_be_fixed_number;
double **DTLZ;
double **ref_points_inside;
double **temp_refpoints;
int *temp_refpoints_pointer;
double **minimum_amount_ref_points;
double **ref_points_normalized;
double *min_ref_points;
double *max_ref_points;
int *index;
double **igb_real_front;
double **igb_algorithm;
double **igb_algorithm_normalized;
double **igb_real_front_normalized;
double *maximum_value;
double *minimum_value;
double *convergence_data;
double *IGD_data;
char *supplied_refpoints_location;
int number_is_feasible;
int number_is_infeasible;
int *feasible_population_sorted_list_index;
int *infeasible_population_sorted_list_index;
int *useless_refpoint_index;
int *usefull_refpoint_index;
int *sort_all_refpoint_index;
int *sort_all_adaptive_refpoint_index;
double **adaptive_refpoints;
double **RRUs;
double **UEs;
int useless_refpoint_number;
int usefull_refpoint_number;
int adaptive_refpoint_number;
int last_gen_adaptive_refpoints_number;
int adaptive_ref_points_inserted;
int adaptive_ref_points_inserted_per_generation;
double *num_div_den;
int IGDfrontsize;
int main (int argc, char **argv)
{
    
    int numberofdivisions;
    int i,j,k;
    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    FILE *fpt5;
    FILE *gp;
    FILE *gp_pc;
    FILE *gp_minus_zmin;
    FILE *gp_normalized;
    FILE *gp_a;
    FILE *gp_dtlz;
    FILE *gp_convergence;
    FILE *gp_IGD;
    FILE *gp_real_front;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;
    population *selection_pop;
    int ref_points_count;
    char ref_points_location;
    
    gp_minus_zmin = popen(GNUPLOT_COMMAND,"w");
    gp_normalized = popen(GNUPLOT_COMMAND,"w");
    gp_a = popen(GNUPLOT_COMMAND,"w");
    gp_real_front = popen(GNUPLOT_COMMAND,"w");
    gp_convergence = popen(GNUPLOT_COMMAND,"w");
    gp_IGD = popen(GNUPLOT_COMMAND,"w");
    if (argc<2)
    {
        printf("\nUsage ./nsga2r random_seed, argc is %d\n",argc);
        exit(1);
    }
    printf("argc is %d\n",argc);
    seed = (double)atof(argv[1]);
    if (seed<=0.0 || seed>=1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    fpt1 = fopen("initial_pop.out","w");
    fpt2 = fopen("final_pop.out","w");
    fpt3 = fopen("best_pop.out","w");
    fpt4 = fopen("all_pop.out","w");
    fpt5 = fopen("params.out","w");
    
    fprintf(fpt1,"# This file contains the data of initial population\n");
    fprintf(fpt2,"# This file contains the data of final population\n");
    fprintf(fpt3,"# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4,"# This file contains the data of all generations\n");
    fprintf(fpt5,"# This file contains information about inputs as read by the program\n");
    printf("\n Enter the problem relevant and algorithm relevant parameters ... ");
    printf("\n Do you want an adaptive NSGA-III simulations? : (2 --> A^2 NSGA-III, 1 --> A NSGA-III, 0 --> NSGA-III) :");
    scanf("%d",&adaptive_nsga);
    printf("\n Enter the population size (a multiple of 4) : ");
    scanf("%d",&popsize);
    if (popsize<4 || (popsize%4)!= 0)
    {
        printf("\n population size read is : %d",popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of generations : ");
    scanf("%d",&ngen);
    if (ngen<1)
    {
        printf("\n number of generations read is : %d",ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of objectives : ");
    scanf("%d",&nobj);
    if (nobj<1)
    {
        printf("\n number of objectives entered is : %d",nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of inequality constraints : ");
    scanf("%d",&ncon);
    if (ncon<0)
    {
        printf("\n number of inequality constraints entered is : %d",ncon);
        printf("\n Wrong number of inequality constraints enetered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of equality constraints : ");
    scanf("%d",&neqcon);
    if (ncon<0)
    {
        printf("\n number of equality constraints entered is : %d",neqcon);
        printf("\n Wrong number of equality constraints enetered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of real variables : ");
    scanf("%d",&nreal);
    if (nreal<0)
    {
        printf("\n number of real variables entered is : %d",nreal);
        printf("\n Wrong number of variables entered, hence exiting \n");
        exit (1);
    }
    if (nreal != 0)
    {
        min_realvar = (double *)malloc(nreal*sizeof(double));
        max_realvar = (double *)malloc(nreal*sizeof(double));
        for (i=0; i<nreal; i++)
        {
            printf ("\n Enter the lower limit of real variable %d : ",i+1);
            scanf ("%lf",&min_realvar[i]);
            printf ("\n Enter the upper limit of real variable %d : ",i+1);
            scanf ("%lf",&max_realvar[i]);
            if (max_realvar[i] < min_realvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
        scanf ("%lf",&pcross_real);
        if (pcross_real<0.0 || pcross_real>1.0)
        {
            printf("\n Probability of crossover entered is : %e",pcross_real);
            printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probablity of mutation of real variables (1/nreal) : ");
        scanf ("%lf",&pmut_real);
        if (pmut_real<0.0 || pmut_real>1.0)
        {
            printf("\n Probability of mutation entered is : %e",pmut_real);
            printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for crossover (5-20): ");
        scanf ("%lf",&eta_c);
        if (eta_c<=0)
        {
            printf("\n The value entered is : %e",eta_c);
            printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for mutation (5-50): ");
        scanf ("%lf",&eta_m);
        if (eta_m<=0)
        {
            printf("\n The value entered is : %e",eta_m);
            printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
            exit (1);
        }
    }
    printf("\n Enter the number of binary variables : ");
    scanf("%d",&nbin);
    if (nbin<0)
    {
        printf ("\n number of binary variables entered is : %d",nbin);
        printf ("\n Wrong number of binary variables entered, hence exiting \n");
        exit(1);
    }
    if (nbin != 0)
    {
        nbits = (int *)malloc(nbin*sizeof(int));
        min_binvar = (double *)malloc(nbin*sizeof(double));
        max_binvar = (double *)malloc(nbin*sizeof(double));
        for (i=0; i<nbin; i++)
        {
            printf ("\n Enter the number of bits for binary variable %d : ",i+1);
            scanf ("%d",&nbits[i]);
            if (nbits[i] < 1)
            {
                printf("\n Wrong number of bits for binary variable entered, hence exiting");
                exit(1);
            }
            printf ("\n Enter the lower limit of binary variable %d : ",i+1);
            scanf ("%lf",&min_binvar[i]);
            printf ("\n Enter the upper limit of binary variable %d : ",i+1);
            scanf ("%lf",&max_binvar[i]);
            if (max_binvar[i] <= min_binvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
        scanf ("%lf",&pcross_bin);
        if (pcross_bin<0.0 || pcross_bin>1.0)
        {
            printf("\n Probability of crossover entered is : %e",pcross_bin);
            printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probability of mutation of binary variables (1/nbits): ");
        scanf ("%lf",&pmut_bin);
        if (pmut_bin<0.0 || pmut_bin>1.0)
        {
            printf("\n Probability of mutation entered is : %e",pmut_bin);
            printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
    }
    if (nreal==0 && nbin==0)
    {
        printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
        exit(1);
    }
    choice=0;
    printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
    scanf("%d",&choice);
    if (choice!=0 && choice!=1)
    {
        printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
        exit(1);
    }
    if (choice==1)
    {
        gp = popen(GNUPLOT_COMMAND,"w");
	gp_pc = popen(GNUPLOT_COMMAND,"w");
        gp_dtlz = popen(GNUPLOT_COMMAND,"w");
	gp_convergence = popen(GNUPLOT_COMMAND,"w");
	gp_IGD = popen(GNUPLOT_COMMAND,"w");

        if (gp==NULL)
        {
            printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
            printf("\n Edit the string to suit your system configuration and rerun the program\n");
            exit(1);
        }
        if (nobj==2)
        {
            printf("\n Enter the objective for X axis display : ");
            scanf("%d",&obj1);
            if (obj1<1 || obj1>nobj)
            {
                printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                exit(1);
            }
            printf("\n Enter the objective for Y axis display : ");
            scanf("%d",&obj2);
            if (obj2<1 || obj2>nobj)
            {
                printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                exit(1);
            }
            obj3 = -1;
        }
        else
        {
            printf("\n #obj > 2, 2D display, 3D display or Parallel Coordinates?, enter 2 for 2D, 3 for 3D and P for Parallel coordinates:");
            scanf("%d",&choice);
            if (choice>25)
            {
                printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
                exit(1);
            }
            if (choice==2)
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&obj1);
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&obj2);
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                obj3 = -1;
            }
            if (choice==3)
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&obj1);
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&obj2);
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                printf("\n Enter the objective for Z axis display : ");
                scanf("%d",&obj3);
                if (obj3<1 || obj3>nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n",obj3);
                    exit(1);
                }
                printf("\n You have chosen 3D display, hence location of eye required \n");
                printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
                scanf("%d",&angle1);
                if (angle1<0 || angle1>180)
                {
                    printf("\n Wrong value for first angle entered, hence exiting \n");
                    exit(1);
                }
                printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
                scanf("%d",&angle2);
                if (angle2<0 || angle2>360)
                {
                    printf("\n Wrong value for second angle entered, hence exiting \n");
                    exit(1);
                }
            }
	    if (choice==4)
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&obj1);
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&obj2);
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                printf("\n Enter the objective for Z axis display : ");
                scanf("%d",&obj3);
                if (obj3<1 || obj3>nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n",obj3);
                    exit(1);
                }
		printf("\n Enter the objective for A axis display : ");
                scanf("%d",&obj4);
                if (obj4<1 || obj4>nobj)
                {
                    printf("\n Wrong value of A objective entered, value entered was %d\n",obj4);
                    exit(1);
                }
                printf("\n You have chosen 4D display \n");

            }
	    if (choice==5)
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&obj1);
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&obj2);
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                printf("\n Enter the objective for Z axis display : ");
                scanf("%d",&obj3);
                if (obj3<1 || obj3>nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n",obj3);
                    exit(1);
                }
		printf("\n Enter the objective for A axis display : ");
                scanf("%d",&obj4);
                if (obj4<1 || obj4>nobj)
                {
                    printf("\n Wrong value of A objective entered, value entered was %d\n",obj4);
                    exit(1);
                }
		printf("\n Enter the objective for B axis display : ");
                scanf("%d",&obj5);
                if (obj5<1 || obj5>nobj)
                {
                    printf("\n Wrong value of A objective entered, value entered was %d\n",obj5);
                    exit(1);
                }
                printf("\n You have chosen 5D display \n");

            }
        }
    }

    printf("\n Input data successfully entered, now performing initialization \n");
    fprintf(fpt5,"\n Population size = %d",popsize);
    fprintf(fpt5,"\n Number of generations = %d",ngen);
    fprintf(fpt5,"\n Number of objective functions = %d",nobj);
    fprintf(fpt5,"\n Number of constraints = %d",ncon);
    fprintf(fpt5,"\n Number of equality constraints = %d",neqcon);
    fprintf(fpt5,"\n Number of real variables = %d",nreal);
    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,min_realvar[i]);
            fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,max_realvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of real variable = %e",pcross_real);
        fprintf(fpt5,"\n Probability of mutation of real variable = %e",pmut_real);
        fprintf(fpt5,"\n Distribution index for crossover = %e",eta_c);
        fprintf(fpt5,"\n Distribution index for mutation = %e",eta_m);
    }
    fprintf(fpt5,"\n Number of binary variables = %d",nbin);
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            fprintf(fpt5,"\n Number of bits for binary variable %d = %d",i+1,nbits[i]);
            fprintf(fpt5,"\n Lower limit of binary variable %d = %e",i+1,min_binvar[i]);
            fprintf(fpt5,"\n Upper limit of binary variable %d = %e",i+1,max_binvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of binary variable = %e",pcross_bin);
        fprintf(fpt5,"\n Probability of mutation of binary variable = %e",pmut_bin);
    }
    fprintf(fpt5,"\n Seed for random number generator = %e",seed);
    bitlength = 0;
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            bitlength += nbits[i];
        }
    }
    fprintf(fpt1,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt2,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt3,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt4,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank\n",nobj,ncon,nreal,bitlength);

    /*Remember numberpointperdim_inside = numberpointperdim-1, if you want to use two layers of reference points*/
    numberpointperdim=(nobj==3)?13:(nobj==4||nobj==5)?7:(nobj==8)?4:(nobj==10)?4:(nobj==15)?3:20;/*numberpointperdim=(nobj==3)?13:(nobj==4)?7:(nobj==5)?7:4;*/
    numberpointperdim_inside=numberpointperdim-1;
    /*popsize=nCn(3-1+12,12)=92(3d),nCn(3-1+6,6)=210(5d),nCn(3-1+8,8)=6435(8d)*/
    numberofdivisions=numberpointperdim-1;
    factorial=fact(numberofdivisions+nobj-1)/(fact(numberofdivisions)*fact(nobj-1));
    factorial_adaptive=fact(nobj)/(fact(1)*fact(nobj-1));
    if (nobj>5)
    	factorial_inside=fact(numberofdivisions+nobj-2)/(fact(numberofdivisions-1)*fact(nobj-1));
    else 	
    	factorial_inside=0;
    last_gen_adaptive_refpoints_number=0;
    elegible_adaptive_ref_points_to_be_fixed_number=0;
    adaptive_ref_points_inserted=0;
    /*printf("factorial %d, factorial_inside %d\n",factorial,factorial_inside);*/
    if (factorial>=popsize)
    {
	printf("The reference points size must be less than the population size. Reduce the number of reference points per dimension! (%d)\n",factorial);
	return (0);
    }
    if (nobj>25)
    {
	printf("NSGA-III algorithm was designed for a number of objectives <=25. Please reduce the dimension!\n");
	return (0);
    }
    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;
    first_front = 0;
  
    a_last_gen=(double *)malloc(nobj*sizeof(double));
    scale_obj_min=(double *)malloc(nobj*sizeof(double));
    scale_obj_max=(double *)malloc(nobj*sizeof(double));
    scale_obj_min_ref=(double *)malloc(nobj*sizeof(int));
    scale_obj_max_ref=(double *)malloc(nobj*sizeof(int));
    a=(double *)malloc(nobj*sizeof(double));
    smin=(double *)malloc(nobj*sizeof(double));
    convergence_data=(double *)malloc(ngen*sizeof(double));
    IGD_data=(double *)malloc(ngen*sizeof(double));
    zmax=(double **)malloc(nobj*sizeof(double*));
    for (i=0;i<nobj;i++)
        zmax[i]=(double *)malloc(nobj*sizeof(double));
    rho=(int *)malloc(nobj*(factorial+factorial_inside)*sizeof(int));
    rho_St=(int *)malloc(nobj*(factorial+factorial_inside)*sizeof(int));
    rho_Fl=(int *)malloc(nobj*(factorial+factorial_inside)*sizeof(int));
    membertoadd=(int *)malloc((2*popsize)*sizeof(int));
    index = (int *)malloc(nobj*sizeof(int));
    memset(index,0,nobj*sizeof(int));
    ref_points_min_rho=(int *)malloc(factorial*sizeof(int));
    ref_points_min_rho_Fl=(int *)malloc(factorial*sizeof(int));
    last_rho_St=(int *)malloc((10*popsize)*sizeof(int));
    last_generation_associated_rho_St=(int *)malloc((10*popsize)*sizeof(int));
    for (i=0;i<nobj*popsize;i++)
    {
	last_rho_St[i]=2147483647;
    }
    index_s=(int *)malloc(nobj*sizeof(int));
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    selection_pop = (population *)malloc(sizeof(population));
    min_ref_points=(double *)malloc(nobj*sizeof(double));
    max_ref_points=(double *)malloc(nobj*sizeof(double));
    useless_refpoint_index = (int *)malloc(popsize*sizeof(int));
    usefull_refpoint_index = (int *)malloc(popsize*sizeof(int));
    sort_all_refpoint_index = (int *)malloc(popsize*sizeof(int));
    sort_all_adaptive_refpoint_index = (int *)malloc(nobj*popsize*sizeof(int));
    num_div_den = (double *)malloc((factorial+factorial_inside)*sizeof(double));

    adaptive_refpoints = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
	adaptive_refpoints[i] = (double *)malloc((factorial+factorial_inside)*sizeof(double));
    }

    minimum_amount_ref_points=(double **)malloc((nobj)*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
		minimum_amount_ref_points[i] = (double *)malloc(sizeof(double));
    }

    temp_refpoints = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
	temp_refpoints[i] = (double *)malloc(((nobj)*(factorial+factorial_inside))*sizeof(double));
    }
    temp_refpoints_pointer=(int *)malloc(((nobj)*(factorial+factorial_inside))*sizeof(int));

    ref_points = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
		ref_points[i] = (double *)malloc(nobj*(factorial+factorial_inside)*sizeof(double));
    }
    adaptive_ref_points_settled = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
		adaptive_ref_points_settled[i] = (double *)malloc(nobj*(factorial+factorial_inside)*sizeof(double));
    }
    adaptive_ref_points_settled_number = (int *)malloc(((nobj)*(factorial+factorial_inside))*sizeof(int));
    DTLZ = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
		DTLZ[i] = (double *)malloc((factorial+factorial_inside)*sizeof(double));
    }

    ref_points_normalized = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
		ref_points_normalized[i] = (double *)malloc((factorial+factorial_inside)*sizeof(double));
    }
    adaptive_refpoints = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
		adaptive_refpoints[i] = (double *)malloc((factorial+factorial_inside)*sizeof(double));
    }
    igb_real_front = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
	igb_real_front[i] = (double *)malloc((popsize)*sizeof(double));
    }
    igb_algorithm = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
	igb_algorithm[i] = (double *)malloc((popsize)*sizeof(double));
    }
    igb_algorithm_normalized = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
	igb_algorithm_normalized[i] = (double *)malloc((popsize)*sizeof(double));
    }
    igb_real_front_normalized = (double **)malloc(nobj*sizeof(double *));
    for (i=0;i<nobj;i++)
    {
	igb_real_front_normalized[i] = (double *)malloc((popsize)*sizeof(double));
    }
    RRUs = (double **)malloc(5*sizeof(double *));
    for (i=0;i<5;i++)
    {
	RRUs[i] = (double *)malloc(6*sizeof(double));
    }
    UEs = (double **)malloc(6*sizeof(double *));
    for (i=0;i<6;i++)
    {
	UEs[i] = (double *)malloc(100*sizeof(double));
    }
    maximum_value = (double *)malloc((nobj)*sizeof(double));
    minimum_value = (double *)malloc((nobj)*sizeof(double));
    w_scalarizing_vector=(double *)malloc(nobj*sizeof(double));
    dist_lf=(int *)malloc((2*popsize)*sizeof(int));
    fronts = (int *)malloc((2*popsize)*sizeof(int));
    feasible_population_sorted_list_index = (int *)malloc((2*popsize)*sizeof(int));
    infeasible_population_sorted_list_index = (int *)malloc((2*popsize)*sizeof(int));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);
    allocate_memory_pop (selection_pop, 2*popsize);

    printf("\n*****************************************************************\nInitialization done, now performing first generation %sNSGA-III\n*****************************************************************\n\n",(adaptive_nsga==2)?"A^2-":(adaptive_nsga==1)?"A-":"");
    printf("Nondominated Sorting Genetic Algorithm version III (%d)\n",(adaptive_nsga==2)?"A^2-":(adaptive_nsga==1)?"A-":"");
    randomize();
    initialize_pop (parent_pop);
    initialize_pop (child_pop);
    initialize_pop (mixed_pop);
    initialize_pop (selection_pop);

    generate_DTLZ1 (nobj-1,1/(double)(numberpointperdim-1));
    display_DTLZ1 ();
    /*onthefly_display_DTLZ1 (gp_dtlz);*/
    if (nobj>5)
    	generate_ref_points_inside (nobj-1,1/(double)(numberpointperdim_inside-1));
    /*dinamic generation of reference points*/
    generate_ref_points (nobj-1,1/(double)(numberpointperdim-1));
    display_refpoints ();
    /*onthefly_display_refpoints (parent_pop, gp_a);
    onthefly_display_one (parent_pop, gp_a);*/
    /*generation of adaptive refpoints cluster*/
    generate_adaptive_ref_points (nobj-1,1.0);

    decode_pop(parent_pop);
    evaluate_pop (parent_pop);
    report_pop (parent_pop, fpt1);
    fprintf(fpt4,"# gen = 1\n");
    report_pop(parent_pop,fpt4);
    printf("\n gen = 1, nobj=%d\n",nobj);
    fflush(stdout);

    /*if (choice!=0)    {
	if (nobj>3)
	{
        	onthefly_display_parallel_coordinates(parent_pop,gp_pc,1);
	}
	else
        	onthefly_display (parent_pop,gp,1,0);
    }*/
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    /*sleep(3);*/
    double temp_IGD=DBL_MAX;
    double IGD_value;
    int temp_gen=0;
    double convergence_value;
    int i1,i2,i3;
    clock_t start;
    for (i=0;i<nobj*(factorial+factorial_inside);i++)
	for (k=0;k<nobj;k++)
    		adaptive_ref_points_settled[k][i]=0;
    for (i=2; i<=ngen; i++)
    {
	clock_t start = clock();
   	 /*if (i==50)
   	 {
		exit(-1);
    	 }*/
        selection (parent_pop, child_pop);/*ok*/
        mutation_pop (child_pop);/*ok*/
        decode_pop(child_pop);/*ok*/
        evaluate_pop(child_pop);/*ok*/
	/*exit(-1);*/
        merge (parent_pop, child_pop, mixed_pop);/*ok*/
        fill_nondominated_sort (selection_pop, mixed_pop, parent_pop,i);
        report_pop(parent_pop,fpt4);
        fflush(fpt4);
        if (choice!=0){
	    if (nobj>3)
	    {
	    	onthefly_display_parallel_coordinates_normalized(parent_pop,gp_pc,i);
	    }
	    else
	    {
	   	onthefly_display (parent_pop,gp,i,1);
	    }
	}
	/*for (i1=0;i1<popsize;i1++)
	{
		for (i2=0;i2<nobj;i2++)
			printf("%e\t",parent_pop->ind[i1].obj[i2]);
		printf("\n");
	}
	exit(-1);*/
	/*if (adaptive_nsga==1 || adaptive_nsga==2)
	{
		printf("Visualization of adaptive reference points\n");
		onthefly_display_refpoints (parent_pop, gp_a);
	}*/
	
	if (dtlz<16)
	{
		IGD_value=IGD(parent_pop);
		if (IGD_value<temp_IGD)
		{
			temp_IGD=IGD_value;
			temp_gen=i;
		}
	}
	convergence_value=convergence_metric();
	if (dtlz<16)
        	printf("\n gen = %d, IGD %e, convergence metric %e\n",i,temp_IGD,convergence_value);
	else
		printf("\n gen = %d",i);
	IGD_data[i]=temp_IGD;
	convergence_data[i]=convergence_value;
	/*sleep(1);*/
    	clock_t end = clock();
    	printf("Runtime %e s\n",(float)(end - start) / CLOCKS_PER_SEC);
	sleep(1);
    }
    /*printf("IGD metric\n");
    for (i=0;i<ngen;i++)
	printf("%e\n",IGD_data[i]);
    printf("convergence metric\n");
    for (i=0;i<ngen;i++)
	printf("%e\n",convergence_data[i]);*/
    if (nobj<=3)
	/*onthefly_display_convergence (gp_convergence,ngen);
	onthefly_display_IGD (gp_IGD,ngen);
    	onthefly_display_real_front (parent_pop,gp_real_front);*/
    if (adaptive_nsga == 1)
    	printf("\nGenerations finished, now reporting solutions (A-NSGA-III)\n");
    else if (adaptive_nsga == 2)
	printf("\nGenerations finished, now reporting solutions (A^2-NSGA-III)\n");
    else
    	printf("\nGenerations finished, now reporting solutions (NSGA-III)\n");
    clock_t end = clock();
    printf("Runtime %e s\n",(float)(end - start) / CLOCKS_PER_SEC);
    report_pop(parent_pop,fpt2);
    report_feasible(parent_pop,fpt3);
    if (dtlz<16)
    {
    	printf("The DTLZ%d Inverted Generational Distance (IGD) for %d dimensions is %e\n",dtlz,nobj,IGD(parent_pop));
    	printf("The best IGD is %e in generation %d\n",temp_IGD,temp_gen);
    }
    if (nreal!=0)
    {
        fprintf(fpt5,"\n Number of crossover of real variable = %d",nrealcross);
        fprintf(fpt5,"\n Number of mutation of real variable = %d",nrealmut);
    }
    if (nbin!=0)
    {
        fprintf(fpt5,"\n Number of crossover of binary variable = %d",nbincross);
        fprintf(fpt5,"\n Number of mutation of binary variable = %d",nbinmut);
    }
    fflush(stdout);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    fclose(fpt4);
    fclose(fpt5);
    if (choice!=0)
    {
        pclose(gp);
        pclose(gp_dtlz);
        pclose(gp_convergence);
        pclose(gp_IGD);
        pclose(gp_pc);

    }
    if (nreal!=0)
    {
        free (min_realvar);
        free (max_realvar);
    }
    if (nbin!=0)
    {
        free (min_binvar);
        free (max_binvar);
        free (nbits);
    }
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2*popsize);
    deallocate_memory_pop (selection_pop, 2*popsize);

    free (parent_pop);
    free (child_pop);
    free (mixed_pop);
    free (selection_pop);
    /*display_fronts ();*/
    printf("\n Routine successfully exited \n");
    return (0);
}
