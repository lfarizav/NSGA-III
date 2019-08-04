/* Test problem definitions */
/* The Copyright belongs to Luis Felipe Ariza Vesga (lfarizav@unal.edu.co). You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for research purposes. All publications which use this code should acknowledge the author. Luis Felipe Ariza Vesga. 
A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems. March, 2019. */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include "global.h"
# include "rand.h"

/* # define sch1 */
/* # define sch2 */
/* # define fon */
/* # define kur */
/* # define pol */
/* # define vnt */
/* # define zdt1 */
/* # define zdt2 */
/* # define zdt3 */
/* # define zdt4*/
 # define c2dtlz2
/*# define dtlz1*/
/* # define zdt5 */
/* # define zdt6  */
/* # define bnh  */
/* # define osy  */
/* # define srn  */
/* # define tnk  */
/* # define ctp1 */
/* # define ctp2 */
/* # define ctp3 */
/* # define ctp4 */
/* # define ctp5 */
/* # define ctp6 */
/* # define ctp7 */
/* # define ctp8 */
/*  # define pipe*/
/*# define pipe4d*/

/*  Test problem crashworthiness
    # of real variables = 5
    # of bin variables = 0
    # of objectives = 3
    # of constraints = 0
    */
#ifdef crashworthiness
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
	obj[0]=1640.2823+2.3573285*xreal[0]+2.3220035*xreal[1]+4.5688768*xreal[2]+7.7213633*xreal[3]+4.4559504*xreal[4];

	obj[1]=6.5856+1.15*xreal[0]-1.0427*xreal[1]+0.9738*xreal[2]
+0.8364*xreal[3]-0.3695*xreal[0]*xreal[3]+0.0861*xreal[0]*xreal[4]+0.3628*xreal[1]*xreal[3]-0.1106*xreal[0]*xreal[0]-0.3437*xreal[2]*xreal[2]+0.1764*xreal[3]*xreal[3];

	obj[2]=-0.0551+0.0181*xreal[0]+0.1024*xreal[1]+0.0421*xreal[2]-0.0073*xreal[0]*xreal[1]+0.024*xreal[1]*xreal[2]-0.0118*xreal[1]*xreal[3]-0.0204*xreal[2]*xreal[3]-0.008*xreal[2]*xreal[4]-0.0241*xreal[1]*xreal[1]+0.0109*xreal[3]*xreal[3];

    return;
}
#endif

/*  Test problem carsideimpact
    # of real variables = 3
    # of bin variables = 0
    # of objectives = 4
    # of constraints = 7
    */
#ifdef waterproblem
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
	dtlz=7;
	obj[0]=106780.37*(xreal[1]+xreal[2])+61704.67;
	obj[1]=3000*xreal[0];
	obj[2]=305700*2289*xreal[1]/pow(0.06*2289,0.65);
	obj[3]=250*2289*exp(-39.75*xreal[1]+9.9*xreal[2]+2.74);
	obj[4]=25*(1.39/(xreal[0]*xreal[1])+4940*xreal[2]-80);
	constr[0]=1-(0.00139/(xreal[0]*xreal[1])+4.94*xreal[2]-0.08);
	constr[1]=1-(0.000306/(xreal[0]*xreal[1])+1.082*xreal[2]-0.0986);
	constr[2]=50000-(12.307/(xreal[0]*xreal[1])+49408.24*xreal[2]+4051.02);
	constr[3]=16000-(2.098/(xreal[0]*xreal[1])+8046.33*xreal[2]-696.71);
	constr[4]=10000-(2.138/(xreal[0]*xreal[1])+7883.39*xreal[2]-705.04);
	constr[5]=2000-(0.417/(xreal[0]*xreal[1])+1721.26*xreal[2]-136.54);
	constr[6]=550-(0.164/(xreal[0]*xreal[1])+631.13*xreal[2]-54.48);
    return;
}
#endif
/*  Test problem waterproblem
    # of real variables = 7
    # of bin variables = 0
    # of objectives = 3
    # of constraints = 10
    */
#ifdef carsideimpact
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
	dtlz=8;
	obj[0] = (1.98+4.9*xreal[0]+6.67*xreal[1]+6.98*xreal[2]+4.01*xreal[3]+1.78*xreal[4]+0.00001*xreal[5]+2.73*xreal[6]);
        obj[1] = (4.72-0.5*xreal[3]-0.19*xreal[1]*xreal[2]);
        obj[2] = (0.5*((10.58-0.674*xreal[0]*xreal[1]-0.67275*xreal[1])+(16.45-0.489*xreal[2]*xreal[6]-0.843*xreal[4]*xreal[5])));

	constr[0] = 1.0-(1.16-0.3717*xreal[1]*xreal[3]-0.0092928*xreal[2]);
	constr[1] = 0.32-(0.261-0.0159*xreal[0]*xreal[1]-0.06486*xreal[0]-0.019*xreal[1]*xreal[6]+0.0144*xreal[2]*xreal[4]+0.0154464*xreal[5]);
	constr[2] = 0.32-(0.214+0.00817*xreal[4]-0.045195*xreal[0]-0.0135168*xreal[0]+0.03099*xreal[1]*xreal[5]-0.018*xreal[1]*xreal[6]+0.007176*xreal[2]+0.023232*xreal[2]-0.00364*xreal[4]*xreal[5]-0.018*xreal[1]*xreal[1]);
	constr[3] = 0.32-(0.74-0.61*xreal[1]-0.031296*xreal[2]-0.031872*xreal[6]+0.227*xreal[1]*xreal[1]);
	constr[4] = 32.0-(28.98+3.818*xreal[2]-4.2*xreal[0]*xreal[1]+1.27296*xreal[5]-2.68065*xreal[6]);
	constr[5] = 32.0-(33.86+2.95*xreal[2]-5.057*xreal[0]*xreal[1]-3.795*xreal[1]-3.795*xreal[1]-3.4431*xreal[6]+1.45728);
	constr[6] = 32.0-(46.36-9.9*xreal[1]-4.4505*xreal[0]);
	constr[7] = 4-(4.72-0.5*xreal[3]-0.19*xreal[1]*xreal[2]);
	constr[8] = 9.9-(10.58-0.674*xreal[0]*xreal[1]-0.67275*xreal[1]);
	constr[9] = 15.7-(16.45-0.489*xreal[2]*xreal[6]-0.843*xreal[4]*xreal[5]);   
    return;
}
#endif
/*  Test problem pipe
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 4
    # of constraints = 0
    */
#ifdef pipe4d
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    if (!normalized){
        obj[0] = cos(xreal[0])*cos(xreal[1]+PI/4);
        obj[1] = cos(xreal[0])*sin(xreal[1]+PI/4);
        obj[2] = sin(xreal[0]);
	obj[3] = -sin(xreal[0]);
    }
    else{
        obj[0] = cos(xreal[0])*cos(xreal[1]+PI/4);
        obj[1] = cos(xreal[0])*sin(xreal[1]+PI/4);
        obj[2] = sin(xreal[0]);
	obj[3] = -sin(xreal[0]);
    }
    return;
}
#endif
#ifdef pipe
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    if (!normalized){
        obj[0] = cos(xreal[0])*cos(xreal[1]+PI/4);
        obj[1] = cos(xreal[0])*sin(xreal[1]+PI/4);
        obj[2] = sin(xreal[0]);
    }
    else{
        obj[0] = cos(xreal[0])*cos(xreal[1]+PI/4);
        obj[1] = cos(xreal[0])*sin(xreal[1]+PI/4);
        obj[2] = sin(xreal[0])/2;
    }
    return;
}
#endif
/*  Test problem SCH1
    # of real variables = 1
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef sch1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    obj[0] = pow(xreal[0],2.0);
    obj[1] = pow((xreal[0]-2.0),2.0);
    return;
}
#endif

/*  Test problem SCH2
    # of real variables = 1
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef sch2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    if (xreal[0]<=1.0)
    {
        obj[0] = -xreal[0];
        obj[1] = pow((xreal[0]-5.0),2.0);
        return;
    }
    if (xreal[0]<=3.0)
    {
        obj[0] = xreal[0]-2.0;
        obj[1] = pow((xreal[0]-5.0),2.0);
        return;
    }
    if (xreal[0]<=4.0)
    {
        obj[0] = 4.0-xreal[0];
        obj[1] = pow((xreal[0]-5.0),2.0);
        return;
    }
    obj[0] = xreal[0]-4.0;
    obj[1] = pow((xreal[0]-5.0),2.0);
    return;
}
#endif

/*  Test problem FON
    # of real variables = n
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef fon
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double s1, s2;
    int i;
    s1 = s2 = 0.0;
    for (i=0; i<nreal; i++)
    {
        s1 += pow((xreal[i]-(1.0/sqrt((double)nreal))),2.0);
        s2 += pow((xreal[i]+(1.0/sqrt((double)nreal))),2.0);
    }
    obj[0] = 1.0 - exp(-s1);
    obj[1] = 1.0 - exp(-s2);
    return;
}
#endif

/*  Test problem KUR
    # of real variables = 3
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef kur
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i;
    double res1, res2;
    res1 = -0.2*sqrt((xreal[0]*xreal[0]) + (xreal[1]*xreal[1]));
    res2 = -0.2*sqrt((xreal[1]*xreal[1]) + (xreal[2]*xreal[2]));
    obj[0] = -10.0*( exp(res1) + exp(res2));
    obj[1] = 0.0;
    for (i=0; i<3; i++)
    {
        obj[1] += pow(fabs(xreal[i]),0.8) + 5.0*sin(pow(xreal[i],3.0));
    }
    return;
}
#endif

/*  Test problem POL
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef pol
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double a1, a2, b1, b2;
    a1 = 0.5*sin(1.0) - 2.0*cos(1.0) + sin(2.0) - 1.5*cos(2.0);
    a2 = 1.5*sin(1.0) - cos(1.0) + 2.0*sin(2.0) - 0.5*cos(2.0);
    b1 = 0.5*sin(xreal[0]) - 2.0*cos(xreal[0]) + sin(xreal[1]) - 1.5*cos(xreal[1]);
    b2 = 1.5*sin(xreal[0]) - cos(xreal[0]) + 2.0*sin(xreal[1]) - 0.5*cos(xreal[1]);
    obj[0] = 1.0 + pow((a1-b1),2.0) + pow((a2-b2),2.0);
    obj[1] = pow((xreal[0]+3.0),2.0) + pow((xreal[1]+1.0),2.0);
    return;
}
#endif

/*  Test problem VNT
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 3
    # of constraints = 0
    */

#ifdef vnt
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
if (nobj==3)
{
    obj[0] = 0.5*(xreal[0]*xreal[0] + xreal[1]*xreal[1]) + sin(xreal[0]*xreal[0] + xreal[1]*xreal[1]);
    obj[1] = (pow((3.0*xreal[0] - 2.0*xreal[1] + 4.0),2.0))/8.0 + (pow((xreal[0]-xreal[1]+1.0),2.0))/27.0 + 15.0;
    obj[2] = 1.0/(xreal[0]*xreal[0] + xreal[1]*xreal[1] + 1.0) - 1.1*exp(-(xreal[0]*xreal[0] + xreal[1]*xreal[1]));
}
if (nobj==4)
{
    obj[0] = 0.5*(xreal[0]*xreal[0] + xreal[1]*xreal[1]) + sin(xreal[0]*xreal[0] + xreal[1]*xreal[1]);
    obj[1] = (pow((3.0*xreal[0] - 2.0*xreal[1] + 4.0),2.0))/8.0 + (pow((xreal[0]-xreal[1]+1.0),2.0))/27.0 + 15.0;
    obj[2] = 1.0/(xreal[0]*xreal[0] + xreal[1]*xreal[1] + 1.0) - 1.1*exp(-(xreal[0]*xreal[0] + xreal[1]*xreal[1]));
    obj[3] = 2.0/(xreal[0]*xreal[0] * xreal[1]*xreal[1] + 1.0) + 1.1*exp(-(xreal[0]*xreal[0] + xreal[1]*xreal[1]));
    /*printf("obj[0] %e, obj[1] %e, obj[2] %e, obj[3] %e\n",obj[0],obj[1],obj[2],obj[3]);*/
}
if (nobj==5)
{
    obj[0] = 0.5*(xreal[0]*xreal[0] + xreal[1]*xreal[1]) + sin(xreal[0]*xreal[0] + xreal[1]*xreal[1]);
    obj[1] = (pow((3.0*xreal[0] - 2.0*xreal[1] + 4.0),2.0))/8.0 + (pow((xreal[0]-xreal[1]+1.0),2.0))/27.0 + 15.0;
    obj[2] = 1.0/(xreal[0]*xreal[0] + xreal[1]*xreal[1] + 1.0) - 1.1*exp(-(xreal[0]*xreal[0] + xreal[1]*xreal[1]));
    obj[3] = 2.0/(xreal[0]*xreal[0] * xreal[1]*xreal[1] + 1.0) + 1.1*exp(-(xreal[0]*xreal[0] + xreal[1]*xreal[1]));
    obj[4] = 1.0/(xreal[0]*xreal[0] / xreal[1]*xreal[1] + 1.0) + 1.1*exp(-(xreal[0]*xreal[0] + xreal[1]*xreal[1]));
    /*printf("obj[0] %e, obj[1] %e, obj[2] %e, obj[3] %e\n",obj[0],obj[1],obj[2],obj[3]);*/
}
    return;
}
#endif

/*  Test problem ZDT1
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double f1, f2, g, h;
    int i;
    f1 = xreal[0];
    g = 0.0;
    for (i=1; i<30; i++)
    {
        g += xreal[i];
    }
    g = 9.0*g/29.0;
    g += 1.0;
    h = 1.0 - sqrt(f1/g);
    f2 = g*h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT2
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double f1, f2, g, h;
    int i;
    f1 = xreal[0];
    g = 0.0;
    for (i=1; i<30; i++)
    {
        g += xreal[i];
    }
    g = 9.0*g/29.0;
    g += 1.0;
    h = 1.0 - pow((f1/g),2.0);
    f2 = g*h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT3
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt3
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double f1, f2, g, h;
    int i;
    f1 = xreal[0];
    g = 0.0;
    for (i=1; i<30; i++)
    {
        g += xreal[i];
    }
    g = 9.0*g/29.0;
    g += 1.0;
    h = 1.0 - sqrt(f1/g) - (f1/g)*sin(10.0*PI*f1);
    f2 = g*h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT4
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt4
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double f1, f2, g, h;
    int i;
    if (nobj==2)
    {
    	f1 = xreal[0];
    	g = 0.0;
    	for (i=1; i<10; i++)
    	{
    	    g += xreal[i]*xreal[i] - 10.0*cos(4.0*PI*xreal[i]);
    	}
    	g += 91.0;
    	h = 1.0 - sqrt(f1/g);
    	f2 = g*h;
    	obj[0] = f1;
    	obj[1] = f2;
    }
    if (nobj==3)
    {
	f1 = xreal[0];
    	g = 0.0;
    	for (i=1; i<10; i++)
    	{
    	    g += xreal[i]*xreal[i] - 10.0*cos(4.0*PI*xreal[i]);
    	}
    	g += 91.0;
    	h = 1.0 - sqrt(f1/g);
    	f2 = g*h;
    	obj[0] = f1;
    	obj[1] = f2;
	obj[2] = f3;
    }
    return;
}
#endif

/*  Test problem ZDT5
    # of real variables = 0
    # of bin variables = 11
    # of bits for binvar1 = 30
    # of bits for binvar2-11 = 5
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt5
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i, j;
    int u[11];
    int v[11];
    double f1, f2, g, h;
    for (i=0; i<11; i++)
    {
        u[i] = 0;
    }
    for (j=0; j<30; j++)
    {
        if (gene[0][j] == 1)
        {
            u[0]++;
        }
    }
    for (i=1; i<11; i++)
    {
        for (j=0; j<4; j++)
        {
            if (gene[i][j] == 1)
            {
                u[i]++;
            }
        }
    }
    f1 = 1.0 + u[0];
    for (i=1; i<11; i++)
    {
        if (u[i] < 5)
        {
            v[i] = 2 + u[i];
        }
        else
        {
            v[i] = 1;
        }
    }
    g = 0;
    for (i=1; i<11; i++)
    {
        g += v[i];
    }
    h = 1.0/f1;
    f2 = g*h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT6
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt6
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double f1, f2, g, h;
    int i;
    f1 = 1.0 - ( exp(-4.0*xreal[0]) ) * pow( (sin(6.0*PI*xreal[0])),6.0 );
    g = 0.0;
    for (i=1; i<10; i++)
    {
        g += xreal[i];
    }
    g = g/9.0;
    g = pow(g,0.25);
    g = 1.0 + 9.0*g;
    h = 1.0 - pow((f1/g),2.0);
    f2 = g*h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem BNH
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef bnh
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    obj[0] = 4.0*(xreal[0]*xreal[0] + xreal[1]*xreal[1]);
    obj[1] = pow((xreal[0]-5.0),2.0) + pow((xreal[1]-5.0),2.0);
    constr[0] = 1.0 - (pow((xreal[0]-5.0),2.0) + xreal[1]*xreal[1])/25.0;
    constr[1] = (pow((xreal[0]-8.0),2.0) + pow((xreal[1]+3.0),2.0))/7.7 - 1.0;
    return;
}
#endif

/*  Test problem OSY
    # of real variables = 6
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 6
    */

#ifdef osy
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    obj[0] = -(25.0*pow((xreal[0]-2.0),2.0) + pow((xreal[1]-2.0),2.0) + pow((xreal[2]-1.0),2.0) + pow((xreal[3]-4.0),2.0) + pow((xreal[4]-1.0),2.0));
    obj[1] = xreal[0]*xreal[0] +  xreal[1]*xreal[1] + xreal[2]*xreal[2] + xreal[3]*xreal[3] + xreal[4]*xreal[4] + xreal[5]*xreal[5];
    constr[0] = (xreal[0]+xreal[1])/2.0 - 1.0;
    constr[1] = 1.0 - (xreal[0]+xreal[1])/6.0;
    constr[2] = 1.0 - xreal[1]/2.0 + xreal[0]/2.0;
    constr[3] = 1.0 - xreal[0]/2.0 + 3.0*xreal[1]/2.0;
    constr[4] = 1.0 - (pow((xreal[2]-3.0),2.0))/4.0 - xreal[3]/4.0;
    constr[5] = (pow((xreal[4]-3.0),2.0))/4.0 + xreal[5]/4.0 - 1.0;
    return;
}
#endif

/*  Test problem SRN
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef srn
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    obj[0] = 2.0 + pow((xreal[0]-2.0),2.0) + pow((xreal[1]-1.0),2.0);
    obj[1] = 9.0*xreal[0] - pow((xreal[1]-1.0),2.0);
    constr[0] = 1.0 - (pow(xreal[0],2.0) + pow(xreal[1],2.0))/225.0;
    constr[1] = 3.0*xreal[1]/10.0 - xreal[0]/10.0 - 1.0;
    return;
}
#endif

/*  Test problem TNK
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef tnk
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    obj[0] = xreal[0];
    obj[1] = xreal[1];
    if (xreal[1] == 0.0)
    {
        constr[0] = -1.0;
    }
    else
    {
        constr[0] = xreal[0]*xreal[0] + xreal[1]*xreal[1] - 0.1*cos(16.0*atan(xreal[0]/xreal[1])) - 1.0;
    }
    constr[1] = 1.0 - 2.0*pow((xreal[0]-0.5),2.0) + 2.0*pow((xreal[1]-0.5),2.0);
    return;
}
#endif

/*  Test problem CTP1
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef ctp1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*exp(-obj[0]/g);
    constr[0] = obj[1]/(0.858*exp(-0.541*obj[0]))-1.0;
    constr[1] = obj[1]/(0.728*exp(-0.295*obj[0]))-1.0;
    return;
}
#endif

/*  Test problem CTP2
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2*PI;
    a = 0.2;
    b = 10.0;
    c = 1.0;
    d = 6.0;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*(1.0  - sqrt(obj[0]/g));
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[0] = exp1/exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP3
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp3
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2*PI;
    a = 0.1;
    b = 10.0;
    c = 1.0;
    d = 0.5;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*(1.0  - sqrt(obj[0]/g));
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[0] = exp1/exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP4
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp4
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2*PI;
    a = 0.75;
    b = 10.0;
    c = 1.0;
    d = 0.5;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*(1.0  - sqrt(obj[0]/g));
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[0] = exp1/exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP5
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp5
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2*PI;
    a = 0.1;
    b = 10.0;
    c = 2.0;
    d = 0.5;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*(1.0  - sqrt(obj[0]/g));
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[0] = exp1/exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP6
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp6
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = 0.1*PI;
    a = 40.0;
    b = 0.5;
    c = 1.0;
    d = 2.0;
    e = -2.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*(1.0  - sqrt(obj[0]/g));
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[0] = exp1/exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP7
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp7
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.05*PI;
    a = 40.0;
    b = 5.0;
    c = 1.0;
    d = 6.0;
    e = 0.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*(1.0  - sqrt(obj[0]/g));
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[0] = exp1/exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP8
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef ctp8
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g*(1.0  - sqrt(obj[0]/g));
    theta = 0.1*PI;
    a = 40.0;
    b = 0.5;
    c = 1.0;
    d = 2.0;
    e = -2.0;
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[0] = exp1/exp2 - 1.0;
    theta = -0.05*PI;
    a = 40.0;
    b = 2.0;
    c = 1.0;
    d = 6.0;
    e = 0.0;
    exp1 = (obj[1]-e)*cos(theta) - obj[0]*sin(theta);
    exp2 = (obj[1]-e)*sin(theta) + obj[0]*cos(theta);
    exp2 = b*PI*pow(exp2,c);
    exp2 = fabs(sin(exp2));
    exp2 = a*pow(exp2,d);
    constr[1] = exp1/exp2 - 1.0;
    return;
}
#endif
/*  Test problem DTLZ1
    # of real variables = 5 + # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=5;
    int k=n_var-nobj+1;
    dtlz=1;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5)-cos(20.0*PI*(xreal[i]-0.5));	
    }
    g = 100*(k+g);
    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g)*0.5;
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=xreal[j];
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= 1 - xreal[aux];
		}
    }
    return;
}
#endif
/*  Test problem inverted DTLZ1
    # of real variables = 5 + # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef idtlz1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=5;
    int k=n_var-nobj+1;
    dtlz=5;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5)-cos(20.0*PI*(xreal[i]-0.5));	
    }
    g = 100*(k+g);			


    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g)*0.5;
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=xreal[j];
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= 1 - xreal[aux];
		}
		obj[i] = (1.0+g)*0.5-obj[i];
    }
    return;
}
#endif
/*  Test problem C1DTLZ1
    # of real variables = 5 + # of objectives -1
    # of bin variables = 0
    # of objectives = n
    # of constraints = 1
    */
#ifdef c1dtlz1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g,temp_constr;
    int n_var=5;
    int k=n_var-nobj+1;
    dtlz=8;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5)-cos(20.0*PI*(xreal[i]-0.5));	
    }
    g = 100*(k+g);
    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g)*0.5;
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=xreal[j];
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= 1 - xreal[aux];
		}
    }
    for (i=0;i<nobj;i++)
    {
	temp_constr=0;
	for (j=0;j<nobj-1;j++)
		temp_constr=temp_constr+obj[j]/0.5;
	constr[i]=1-obj[nobj-1]/0.6-temp_constr;
    }
    /*for (i=0;i<nobj-2;i++)
    {
	temp_constr+=obj[i]/0.5;
    }
    constr[0]=1-obj[nobj-1]/0.6-temp_constr;*/
    /*printf("temp_constr is %e\n",constr[0]);*/
    return;
}
#endif
/*  Test problem Inverted DTLZ2
    # of real variables = 10 + # of objectives -1 = 12-13-14
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef idtlz2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    int n_var=10;
    int k=n_var-nobj+1;
    double g;
    dtlz=6;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5);	
    }
    for (i=0;i<nobj;i++)
    {
	obj[i]=1.0+g;	
    }   
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(xreal[j]*0.5*PI);	
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(xreal[aux]*0.5*PI);
		}	
    } 
    for (i=0;i<nobj;i++)
    {
	
	if (i<nobj-1)
		obj[i] = ((1.0+g)-(obj[i]*obj[i]*obj[i]*obj[i]));
	else
		obj[i] = ((1.0+g)-(obj[i]*obj[i]));
    } 
    return;
}
#endif
/*  Test problem DTLZ2
    # of real variables = 10 + # of objectives -1 = 12-13-14
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    int n_var=10;
    int k=n_var-nobj+1;
    double g;
    dtlz=2;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5);	
    }
    for (i=0;i<nobj;i++)
    {
	obj[i]=1+g;	
    }   
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(xreal[j]*0.5*PI);
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(xreal[aux]*0.5*PI);
		}
    } 
    return;
}
#endif
/*  Test problem C2DTLZ2
    # of real variables = 10 + # of objectives -1 = 12-13-14
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 1
    */
#ifdef c2dtlz2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    int n_var=10;
    int k=n_var-nobj+1;
    double g;
    double temp_constr1=0;
    double temp_constr2=0;
    double temp_min1=0;
    double r;
    dtlz=10;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5);	
    }
    for (i=0;i<nobj;i++)
    {
	obj[i]=1+g;	
    }   
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(xreal[j]*0.5*PI);
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(xreal[aux]*0.5*PI);
		}
    } 
    /*Follows the constraint*/
    if (nobj==3)
	r=0.4;
    else
	r=0.5;
    temp_min1=DBL_MAX;
    for (i=0;i<nobj;i++)
    {
	temp_constr1=0;
	temp_constr1=temp_constr1+pow(obj[i]-1,2.0);
	for (j=0;j<nobj;j++)
	{
		if (i!=j)
			temp_constr1=temp_constr1+pow(obj[j],2.0);
	}
	temp_constr1=temp_constr1-pow(r,2.0);
	if (temp_constr1<temp_min1)
		temp_min1=temp_constr1;
    }    

    temp_constr2=0;
    for (i=0;i<nobj;i++)
    {
	temp_constr2=temp_constr2+pow(obj[i]-1/sqrt(nobj),2.0);
    }
    temp_constr2=temp_constr2-pow(r,2.0);

    constr[0]=(temp_min1<temp_constr2)?-temp_min1:-temp_constr2;
    return;
}
#endif
/*  Test problem DTLZ3
    # of real variables = 10 + # of objectives -1 = 12-13-14
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz3
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=10;
    int k=n_var-nobj+1;
    dtlz=3;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5)-cos(20.0*PI*(xreal[i]-0.5));	
    }
    g = 100*(k+g);
    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g);
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(xreal[j]*0.5*PI);
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(xreal[aux]*0.5*PI);
		}
    }
    return;
}
#endif
/*  Test problem C1DTLZ3
    # of real variables = 10 + # of objectives -1 = 12-13-14
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef c1dtlz3
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    double temp_constr=0;
    double r[15]={0,0,9,0,12.5,0,0,12.5,0,15,0,0,0,0,15};
    /*double M[5]={3,5,8,10,15};*/
    int n_var=10;
    int k=n_var-nobj+1;
    dtlz=9;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5)-cos(20.0*PI*(xreal[i]-0.5));	
    }
    g = 100*(k+g);
    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g);
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(xreal[j]*0.5*PI);
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(xreal[aux]*0.5*PI);
		}
    }
    /*for (i=0;i<nobj;i++)
    {*/
	temp_constr=0;
	for (j=0;j<nobj;j++)
	{
		temp_constr=temp_constr+pow(obj[j],2);
	}
    temp_constr=temp_constr-16;
    temp_constr=temp_constr-pow(r[nobj-1],2);
    constr[0]=temp_constr*temp_constr;
    /*}*/
    /*for (i=0;i<nobj;i++)
    {
	temp_constr=0;
	for (j=0;j<nobj-1;j++)
		temp_constr=temp_constr+obj[j]/0.5;
	constr[i]=1-obj[nobj-1]/0.6-temp_constr;
    }*/
    return;
}
#endif
/*  Test problem DTLZ4
    # of real variables = 10+ # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz4
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=10;
    int k=n_var-nobj+1;
    int alpha=100;
    dtlz=4;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5);	
    }
    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g);
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(pow(xreal[j],alpha)*0.5*PI);
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(pow(xreal[aux],alpha)*0.5*PI);
		}
    }
    return;
}
#endif
/*  Test problem DTLZ5
    # of real variables = 10+ # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz5
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=10;
    int k=n_var-nobj+1;
    double tetha[nobj-1];
    dtlz=5;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=(xreal[i]-0.5)*(xreal[i]-0.5);	
    }
    double t = PI/(4.0*(1+g));
    tetha[0]=xreal[0]*0.5*PI;
    for (i=1;i<nobj-1;i++)
    {
	tetha[i]=t*(1.0+2.0*g*xreal[i]);
    }
    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g);
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(tetha[j]);
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(tetha[aux]);
		}
    }
    return;
}
#endif
/*  Test problem DTLZ6
    # of real variables = 10+ # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz6
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=10;
    int k=n_var-nobj+1;
    int alpha=100;
    double tetha[nobj-1];
    dtlz=6;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=pow(xreal[i],0.1);	
    }
    int t = PI/(4.0*(1+g));
    tetha[0]=xreal[0]*0.5*PI;
    for (i=1;i<nobj-1;i++)
    {
	tetha[i]=t*(2.0*g*xreal[i]);
    }
    for (i=0;i<nobj;i++)
    {
	obj[i]=(1.0+g);
    }
    for (i=0;i<nobj;i++)
    {
	for (j=0;j<nobj-(i+1);j++)
		obj[i]*=cos(tetha[j]);
		if (i!=0)
		{
			aux = nobj-(i+1);
			obj[i] *= sin(tetha[aux]);
		}
    }
    return;
}
#endif
/*  Test problem DTLZ7
    # of real variables = 20+ # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz7
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=20;
    int k=n_var-nobj+1;
    dtlz=7;
    g=0.0;
    for (i=n_var-k;i<n_var;i++)
    {
	g+=xreal[i];	
    }
    g=1.0+(9.0*g)/k;
    for (i=0;i<nobj-1;i++)
    {
	obj[i]=xreal[i];	
    }
    double h=0;


    for (i=0;i<nobj-1;i++)
    {
	h+=(obj[i]/(1.0+g))*(1.0+sin(3.0*PI*obj[i]));
    }
    h=nobj-h;
    obj[nobj-1]=(1.0+g)*h;

    return;
}
#endif
/*  Test problem DTLZ7_g
    # of real variables = 20+ # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
#ifdef dtlz7_g
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized)
{
    int i,j,aux;
    double g;
    int n_var=20;
    int k=n_var-nobj+1;
    dtlz=7;
    g=1.0;
    for (i=0;i<nobj-1;i++)
    {
	obj[i]=xreal[i];	
    }
    double h=0;
    for (i=0;i<nobj-1;i++)
    {
	h+=(obj[i]/(1.0+g))*(1.0+sin(3.0*PI*obj[i]));
    }
    h=nobj-h;
    obj[nobj-1]=(1.0+g)*h;

    return;
}
#endif
