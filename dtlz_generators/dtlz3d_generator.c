#include <stdio.h>
#include <math.h>
# define PI 3.14159265358979

int main()
{
int i,j,k,temp=0;
double x1,x2,x3;

	for (x1=0,i=0; x1<=0.26; x1+=0.065,i++)
	{
		for (x2=0,j=0; x2<=0.26; x2+=0.065,j++)
		{
			x3 = 2*(3 - 0.5*(x1*(1+sin(3*PI*x1)) + x2*(1+sin(3*PI*x2))));
			printf("%e\t%e\t%e\n",temp,x1,x2,x3);
			temp++;
		}
		for (x2=0.64,j=0; x2<0.88; x2+=0.056,j++)
		{
			x3 = 2*(3 - 0.5*(x1*(1+sin(3*PI*x1)) + x2*(1+sin(3*PI*x2))));
			printf("%e\t%e\t%e\n",temp,x1,x2,x3);
			temp++;
		}
	}
	for (x1=0.64,i=0;x1<0.88; x1+=0.056,i++)
	{
		for (x2=0,j=0; x2<=0.26; x2+=0.065,j++)
		{
			x3 = 2*(3 - 0.5*(x1*(1+sin(3*PI*x1)) + x2*(1+sin(3*PI*x2))));
			printf("%e\t%e\t%e\n",temp,x1,x2,x3);
			temp++;
		}
		for (x2=0.64,j=0; x2<0.88; x2+=0.056,j++)
		{
			x3 = 2*(3 - 0.5*(x1*(1+sin(3*PI*x1)) + x2*(1+sin(3*PI*x2))));
			printf("%e\t%e\t%e\n",temp,x1,x2,x3);
			temp++;
		}
	}
return (0);
}
