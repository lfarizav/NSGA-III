/* 
 *@Coded by Luis Felipe Ariza Vesga
 *@Universidad Nacional de Colombia
 *@EURECOM
 *@2017
 *@lfarizav@unal.edu.co, ariza@eurecom.fr
printf(".-. .-. .-. .-.   .--.   .-.                
printf("| } { | |  \{ |  / {} \  } |                        
printf("\ `-' / | }\  { /  /\  \ } '--.                     
printf(" `---'  `-' `-' `-'  `-' `----'                     
printf(",----, ,-, ,-. ,---.  .----. .----.  .---.  .-.  .-.    
printf("} |__} | } { | } }}_} } |__} | }`-' / {-. \ }  \/  {    
printf("} '__} \ `-' / | } \  } '__} | },-. \ '-} / | {  } |    
printf("`----'  `---'  `-'-'  `----' `----'  `---'  `-'  `-' 
 */

/* Include standard headers: */
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include<stdio.h>
#include<string.h>
#include<pthread.h>
#include<stdlib.h>
#include<unistd.h>


#define DEFAULTGENELENGTH 63
#define FITNESS 1
# define GNUPLOT_COMMAND "gnuplot -persist"


int main (){
FILE *gp;
gp = popen(GNUPLOT_COMMAND,"w");
FILE *fpt;
fpt = fopen("plot.out","w");
double step;
int i,j,k,l,number, nobj;
number=0;
step=0.5;
    nobj=5;
    
    if (nobj==3){
		for (k=0;k<=(int)(1/step);k++){
			for (l=0;l<=(int)(1/step)-k;l++){
				number++;
				printf("x %e\t y %e\t z %e\n",step*k,step*l,(1-step*(k+l)));
				fprintf(fpt,"%e\t%e\t%e\n",step*k,step*l,(1-step*(k+l)));
			
			}
		}
		printf("number is %d\n",number);
		fprintf(gp,"set title 'Generation #'\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n");
    }
    if (nobj==4){
        for (j=0;j<=(int)(1/step);j++){
            for (k=0;k<=(int)(1/step)-j;k++){
                for (l=0;l<=(int)(1/step)-k-j;l++){
                    number++;
                    printf("x %e\t y %e\t z %e\t a %e \n",step*j,step*k,step*l,(1-step*(j+k+l)));
                    fprintf(fpt,"%e\t%e\t%e\t%e\n",step*j,step*k,step*l,(1-step*(j+k+l)));
                
                }
            }
        }
        printf("number is %d\n",number);
        fprintf(gp,"set title 'Generation #'\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n");
    }
    if (nobj==5){
        for(i=0;i<=(int)(1/step);i++){
            for (j=0;j<=(int)(1/step)-i;j++){
                for (k=0;k<=(int)(1/step)-j-i;k++){
                    for (l=0;l<=(int)(1/step)-k-j-i;l++){
                        number++;
                        printf("x %e\t y %e\t z %e\t a %e\t b %e \n",step*i,step*j,step*k,step*l,(1-step*(i+j+k+l)));
                        fprintf(fpt,"%e\t%e\t%e\t%e\t%e\n",step*i,step*j,step*k,step*l,(1-step*(i+j+k+l)));
                    
                    }
                }
            }
        }
        printf("number is %d\n",number);
        fprintf(gp,"set title 'Generation #'\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n");
    }
}
