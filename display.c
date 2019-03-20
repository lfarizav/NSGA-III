/* Routines to display the population information using gnuplot */
/* The Copyright belongs to Luis Felipe Ariza Vesga (lfarizav@unal.edu.co). You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for research purposes. All publications which use this code should acknowledge the author. Luis Felipe Ariza Vesga. 
A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems. March, 2019. */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <unistd.h>
# include <float.h>

# include "global.h"
# include "rand.h"

/* Function to display the current population for the subsequent generation */
void onthefly_display (population *pop, FILE *gp, int ii, int normalization)
{
    int i,j;
    int flag;
    FILE *fpt;
    fpt = fopen("plot.out","w");
    flag = 0;

    for (i=0; i<popsize; i++)
    {
            if (choice!=3)
                fprintf(fpt,"%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1]);
            else
	    {
			fprintf(fpt,"%e\t%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1],pop->ind[i].obj[obj3-1]);
	    }
            fflush(fpt);
            flag = 1;
        /*}*/
    }
    if (flag==0)
    {
        printf("\n No feasible solutions in this pop, hence no display");
    }
    else
    {
        if (choice!=3)
            fprintf(gp,"set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n",ii);
        else
            fprintf(gp,"set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n",ii,angle1,angle2);
        fflush(gp);
    }
    fclose(fpt);
    return;
}
/* Function to display the real front */
void onthefly_display_real_front (population *pop,FILE *gp)
{
    int i,j;
    FILE *fpt;
    fpt = fopen("plot.out","w");

    for (i=0; i<popsize; i++)
    {
	fprintf(fpt,"%e\t%e\t%e\n",(pop->ind[i].obj[obj1-1]),(pop->ind[i].obj[obj2-1]),(pop->ind[i].obj[obj3-1]));
    }
    if (dtlz==1)
    fprintf(gp,"set title 'Real front'\n set view %d,%d\nunset key\nsplot '/home/lfarizav/Dropbox/EURECOM/nsga2-gnuplot-v1.1.6/real_front/DTLZ1-3-PF.txt' w points pointtype 6 pointsize 1, 'plot.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    if (dtlz==2)
    fprintf(gp,"set title 'Real front'\n set view %d,%d\nunset key\nsplot '/home/lfarizav/Dropbox/EURECOM/nsga2-gnuplot-v1.1.6/real_front/DTLZ2-3-PF.txt' w points pointtype 6 pointsize 1, 'plot.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    if (dtlz==3)
    fprintf(gp,"set title 'Real front'\n set view %d,%d\nunset key\nsplot '/home/lfarizav/Dropbox/EURECOM/nsga2-gnuplot-v1.1.6/real_front/DTLZ3-3-PF.txt' w points pointtype 6 pointsize 1, 'plot.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    if (dtlz==4)
    fprintf(gp,"set title 'Real front'\n set view %d,%d\nunset key\nsplot '/home/lfarizav/Dropbox/EURECOM/nsga2-gnuplot-v1.1.6/real_front/DTLZ4-3-PF.txt' w points pointtype 6 pointsize 1, 'plot.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    if (dtlz==5)
    fprintf(gp,"set title 'Real front'\n set view %d,%d\nunset key\nsplot '/home/lfarizav/Dropbox/EURECOM/nsga2-gnuplot-v1.1.6/real_front/IDTLZ1-3-PF.txt' w points pointtype 6 pointsize 1, 'plot.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    fflush(gp);
    fflush(fpt);
    fclose(fpt);
    return;
}
/* Function to display the current population after substract zmin */
void onthefly_display_minus_zmin (population *pop, FILE *gp)
{
    int i;
    FILE *fpt;
    fpt = fopen("plot_minus_zmin.out","w");
    for (i=0; i<popsize; i++)
    {
            fprintf(fpt,"%e\t%e\t%e\n",pop->ind[i].obj_minus_zmin[0],pop->ind[i].obj_minus_zmin[1],pop->ind[i].obj_minus_zmin[2]);
            fflush(fpt);
    }
            fprintf(gp,"set title 'Generation minus zmin'\n set view %d,%d\n unset key\n splot 'plot_minus_zmin.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
        
    fflush(gp);
    fclose(fpt);

    onthefly_display_zmax (pop,gp);
    return;
}

/* Function to display the current population normalized */
void onthefly_display_normalized (population *pop, FILE *gp, int archieve_and_front_sizes)
{
    int i;
    FILE *fpt;
    fpt = fopen("plot_normalized.out","w");

    for (i=0; i<archieve_and_front_sizes; i++)
    {
        if (pop->ind[i].constr_violation==0)
        {
                fprintf(fpt,"%e\t%e\t%e\n",pop->ind[i].obj_normalized[obj1-1],pop->ind[i].obj_normalized[obj2-1],pop->ind[i].obj_normalized[obj3-1]);
            fflush(fpt);
        }
    }
            fprintf(gp,"set title 'Generation normalized'\n set view %d,%d\n unset key\n splot 'plot_normalized.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    fflush(gp);
    fclose(fpt);
    return;
}
/* Function to display the a matrix */
void onthefly_display_a (population *pop, FILE *gp)
{
    int i;
    int flag;
    FILE *fpt;
    fpt = fopen("plot_a.out","w");

    fprintf(fpt,"%e\t%d\t%d\n",a[0],0,0);
    fprintf(fpt,"%d\t%e\t%d\n",0,a[1],0);
    fprintf(fpt,"%d\t%d\t%e\n",0,0,a[2]);
    fprintf(fpt,"%e\t%d\t%d\n",a[0],0,0);
    fflush(fpt);
    flag = 1;
    fprintf(gp,"set view %d,%d\n unset key\n replot 'plot_a.out' w lines pointtype 6 pointsize 1\n",angle1,angle2);
    fflush(gp);
    fclose(fpt);
    return;
}
/* Function to display the plane x+y+z=1 */
void onthefly_display_one (population *pop, FILE *gp)
{
    int i,j;
    int flag;
    FILE *fpt;
    fpt = fopen("plot_one.out","w");

    fprintf(fpt,"%d\t%d\t%d\n",1,0,0);
    fprintf(fpt,"%d\t%d\t%d\n",0,1,0);
    fprintf(fpt,"%d\t%d\t%d\n",0,0,1);
    fprintf(fpt,"%d\t%d\t%d\n",1,0,0);

    fflush(fpt);
    flag = 1;
    fprintf(gp,"replot 'plot_one.out' w lines pointtype 6 pointsize 1\n");
    fflush(gp);
    fclose(fpt);
    return;
}
/* Function to display the reference points of the inside layer */
void onthefly_display_inside (population *pop, FILE *gp)
{
    int i,j;
    int flag;
    FILE *fpt;
    fpt = fopen("plot_inside.out","w");

    for (i=factorial; i<factorial+factorial_inside; i++)
    {
	for (j=0;j<nobj;j++)
	{
        	fprintf(fpt,"%e\t",ref_points[j][i]);
	}
	fprintf(fpt,"\n");
    }
    fprintf(fpt,"%e\t%e\t%e\t",ref_points[j][i]);
    fflush(fpt);
    flag = 1;
    fprintf(gp,"set view %d,%d\n unset key\n replot 'plot_inside.out' w lines pointtype 6 pointsize 1\n",angle1,angle2);
    fflush(gp);
    fclose(fpt);
    return;
}
/* Function to display the reference points */
void onthefly_display_refpoints (population *pop, FILE *gp)
{
    int i,j;
    int flag;
    FILE *fpt;
    fpt = fopen("plot_rf.out","w");
    for (i=0; i<factorial+factorial_inside; i++)
    {
	for (j=0;j<nobj;j++)
	{
        	fprintf(fpt,"%e\t",ref_points[j][i]);
	}
	fprintf(fpt,"\n");
    }
    fflush(fpt);

    if (adaptive_nsga==1)
    {
 	   fprintf(gp,"set title 'Adaptive Reference Points'\n set view %d,%d\n unset key\n splot 'plot_rf.out' w points pointtype 6 pointsize 1\n",angle1,angle2);

    }
    if (adaptive_nsga==2)
 	   fprintf(gp,"set title 'Effective Adaptive Reference Points'\n set view %d,%d\n unset key\n splot 'plot_rf.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    if (adaptive_nsga==0)
 	   fprintf(gp,"set title 'Reference Points (Das and Dennis approach)'\n set view %d,%d\n unset key\n splot 'plot_rf.out' w points pointtype 6 pointsize 1\n",angle1,angle2);

    printf("adaptive_nsga = %d\n",adaptive_nsga);
    fflush(gp);
    fclose(fpt);
    return;
}
/* Function to display the true pareto of the DTLZ1 problem */
void onthefly_display_DTLZ1 (FILE *gp_dtlz)
{
    int i,j;
    int flag;
    FILE *fpt;
    fpt = fopen("plot_dtlz.out","w");
    for (i=0; i<factorial+factorial_inside; i++)
    {
	for (j=0;j<nobj;j++)
	{
        	fprintf(fpt,"%e\t",DTLZ[j][i]);
	}
	fprintf(fpt,"\n");
    }
    fflush(fpt);
    if (adaptive_nsga==1)
    {
 	   fprintf(gp_dtlz,"set title 'DTLZ1'\n set view %d,%d\n unset key\n splot 'plot_dtlz.out' w points pointtype 6 pointsize 1\n",angle1,angle2);

    }
    if (adaptive_nsga==2)
 	   fprintf(gp_dtlz,"set title 'DTLZ2'\n set view %d,%d\n unset key\n splot 'plot_dtlz.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    if (adaptive_nsga==0)
 	   fprintf(gp_dtlz,"set title 'DTLZ0'\n set view %d,%d\n unset key\n splot 'plot_dtlz.out' w points pointtype 6 pointsize 1\n",angle1,angle2);
    fflush(gp_dtlz);
    fclose(fpt);
    return;
}

/* Function to display reference lines */
void onthefly_display_reflines (population *pop, FILE *gp)
{
    int i,j,k;
    int flag;
    FILE *fpt;
    fpt = fopen("plot_lines.out","w");
    for (i=0; i<factorial+factorial_inside; i++)
    {
	for (k=0;k<nobj;k++)
	{
	fprintf(fpt,"%d\t",0);
	}
	fprintf(fpt,"\n");
	for (k=0;k<nobj;k++)
	{
        fprintf(fpt,"%e\t",ref_points[k][i]);
	}
	fprintf(fpt,"\n");
    }
    fflush(fpt);
    flag = 1;

        fprintf(gp,"set view %d,%d\n unset key\n replot 'plot_lines.out' w lines pointtype 6 pointsize 1\n",angle1,angle2);
	
    
    fflush(gp);
    fclose(fpt);
    return;
}
/* Function to display the current population for the subsequent generation */
void onthefly_display_association (population *pop, FILE *gp)
{
    int i,j,k;
    FILE *fpt;
    fpt = fopen("plot_lines_association.out","w");
    for (i=0; i<popsize; i++)
    {
	for (k=0;k<nobj;k++)
	{
	fprintf(fpt,"%e\t",0.0);
	}
	fprintf(fpt,"\n");

	for (k=0;k<nobj;k++)
	{
	printf("i %d\t ref %d\t rp %e\t w %e\n",i,pop->ind[i].associatedref,ref_points[k][pop->ind[i].associatedref],pop->ind[i].w);
        fprintf(fpt,"%e\t",ref_points[k][pop->ind[i].associatedref]*pop->ind[i].w);
	}
	fprintf(fpt,"\n");
	printf("\n\n");

	for (k=0;k<nobj;k++)
	{
        fprintf(fpt,"%e\t",pop->ind[i].obj_normalized[k]);
	}
	fprintf(fpt,"\n");

	for (k=0;k<nobj;k++)
	{
        fprintf(fpt,"%e\t",ref_points[k][pop->ind[i].associatedref]*pop->ind[i].w);
	}
	fprintf(fpt,"\n");

    }
    fflush(fpt);

        fprintf(gp,"set view %d,%d\n unset key\n replot 'plot_lines_association.out' w lines pointtype 6 pointsize 1\n",angle1,angle2);
	
    
    fflush(gp);
    fclose(fpt);
    return;
}
/* Function to display the zmax matrix */
void onthefly_display_zmax (population *pop, FILE *gp)
{
    int i;
    FILE *fpt;
    fpt = fopen("plot_zmax.out","w");
    for (i=0; i<nobj; i++)
    {
        fprintf(fpt,"%e\t%e\t%e\n",zmax[0][i],zmax[1][i],zmax[2][i]);
        fflush(fpt);
    }
    fprintf(fpt,"%e\t%e\t%e\n",zmax[0][0],zmax[1][0],zmax[2][0]);/*used to close the triangle*/
    fflush(fpt);
    fprintf(gp,"set view %d,%d\n unset key\n replot 'plot_zmax.out' w lines pointtype 6 pointsize 1\n",angle1,angle2);
    fflush(gp);
    fclose(fpt);
    return;
}

/* Function to display the current population in parallel coordinates for the subsequent generation in more than 3 dimmensions	*/
void onthefly_display_parallel_coordinates (population *pop, FILE *gp_pc, int ii)
{
    int i,j,k;
    for (i=0;i<nobj;i++)
    {
	scale_obj_min[i]=DBL_MAX;
    }
    memset(scale_obj_max,0,nobj*sizeof(double));
    FILE *fpt, *fpt_pc;
    fpt_pc = fopen("plot_pc.out","w");
    fpt = fopen("plot.out","w");
    for (j=0; j<popsize; j++)
    {
       	find_min_from_functions(&(pop->ind[j]),j,1);
        find_max_from_functions(&(pop->ind[j]),j,1);
    }
    for (i=0; i<popsize; i++)
    {
		for (j=0;j<nobj;j++)
		{
                	fprintf(fpt_pc,"%d\t%e\n",j+1,(pop->ind[i].obj[j]-scale_obj_min[j])/(scale_obj_max[j]-scale_obj_min[j]));
		}
		for (j=0;j<=nobj;j++)
		{
			fprintf(fpt_pc,"\n");
		}
                fflush(fpt_pc);
    }

    for (i=0; i<popsize; i++)
    {
	    for (j=0;j<nobj;j++)
	    {
			fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
	    }
	    fprintf(fpt,"\n");
            fflush(fpt);
    }
    fprintf(gp_pc,"set title 'Parallel Coordinates'\n unset key\n unset xtics\n set tmargin 6\n");
    fprintf(gp_pc,"set x2tics rotate by 90 offset 0 mirror ( ");
    for (k=1;k<nobj;k++)
    {
			fprintf(gp_pc,"'objective %d' %d,",k,k);
    }
    fprintf(gp_pc,"'objective %d' %d)\n plot for [i=0:%d] 'plot_pc.out' i i u 1:2 w linesp notitle\n",k,k,popsize);

    fflush(gp_pc);
    fclose(fpt_pc);
    fclose(fpt);
    return;
}
void display_pop (population *pop)
{
    int i;
    printf("The popsize is %d\n",popsize);
    for (i=0; i<popsize; i++)
    {
        display_pop_ind_xreal (&(pop->ind[i]),i);
    }
    for (i=0; i<popsize; i++)
    {
        display_pop_ind_obj (&(pop->ind[i]),i);
    }
    for (i=0; i<popsize; i++)
    {
        display_pop_ind_obj_minus_zmin (&(pop->ind[i]),i);
    }
    for (i=0; i<popsize; i++)
    {
        display_pop_ind_obj_normalized (&(pop->ind[i]),i);
    }
}
void display_pop_ind_obj (individual *ind, int popsizeline)
{
    int i;
    double temp=0;
    printf("%d\t",popsizeline);
    for (i=0;i<nobj;i++)
    {
        printf("%e\t",ind->obj[i]);
	temp=temp+ind->obj[i];	
    }
    printf("\n");
    return;
}
void display_pop_ind_obj_minus_zmin (individual *ind, int popsizeline)
{
 	int i;
	printf("%d\t",popsizeline);
 	for (i=0;i<nobj;i++)
	{
        	printf("%e\t",ind->obj_minus_zmin[i]);
	}
	printf("\n");
	return;
}
void display_pop_ind_obj_normalized (individual *ind, int popsizeline)
{
 	int i;
	printf("%d\t",popsizeline);
 	for (i=0;i<nobj;i++)
	{
		printf("%e\t",ind->obj_normalized[i]);
	}
	printf("\n");
        return;
}
void display_pop_ind_xreal (individual *ind, int popsizeline)
{
    int i;
	printf("%d\t",popsizeline);
        for (i=0;i<nreal;i++)
        {
            printf("%e\t",ind->xreal[i]);
        }
        printf("\n");
	return;
}
void display_refpoints ()
{
    int i,j;
    printf("Visualization of reference points:\n");
    printf("factorial: %d\n",factorial);
    if (nobj>1)
    {
    printf("factorial_inside: %d\n",factorial_inside);
    printf("last_gen_adaptive_refpoints_number: %d\n",last_gen_adaptive_refpoints_number);
    printf("adaptive_refpoint_number: %d\n",adaptive_refpoint_number);
    }
    for (i=0; i<factorial+factorial_inside+adaptive_refpoint_number; i++)
    {
	printf("%d\t",i);
	for (j=0;j<nobj;j++)
	{
		printf("%e\t",ref_points[j][i]);
	}
	printf("\n");
    }
}
void display_DTLZ1 ()
{
    int i,j;
    printf("Visualization of DTLZ:\n");

    for (i=0; i<factorial+factorial_inside; i++)
    {
	for (j=0;j<nobj;j++)
	{
		printf("%e\t",DTLZ[j][i]);
	}
	printf("\n");
    }
}
void display_refpoints_normalized ()
{
    int i,j;
    printf("Visualization of reference points normalized:\n");
    printf("factorial: %d\n",factorial);
    if (nobj>5)
    {
    	printf("factorial_inside: %d\n",factorial_inside);
    	printf("factorial_adaptive: %d\n",factorial_adaptive);
    	printf("last_gen_adaptive_refpoints_number: %d\n",last_gen_adaptive_refpoints_number);
    }
    for (i=0; i<factorial+factorial_inside; i++)
    {
	printf("%d\t",i);
	for (j=0;j<nobj;j++)
	{
		printf("%e\t",ref_points_normalized[j][i]);
	}
	printf("\n");
    }
}
void display_fronts ()
{
    int i,j;
    printf("Visualization of algorithm fronts\n");
    for (i=0; i<popsize; i++)
    {
	printf("%d\t",i);
	for (j=0;j<nobj;j++)
	{
		printf("%e\t",igb_algorithm[j][i]);
	}
	printf("\n");
    }
    printf("0. visualization of real fronts \n");
    for (i=0; i<popsize; i++)
    {
	printf("%d\t",i);
	for (j=0;j<nobj;j++)
	{
		printf("%e\t",igb_real_front[j][i]);
	}
	printf("\n");
    }
}
