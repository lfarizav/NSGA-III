# NSGA-III
Non official NSGA-III

This project is based on the NSGA-II official C code from Kanpur Genetic Algorithms Laboratory. It is an extention of NSGA-II to obtain the NSGA-III algorithm. The difference between those algorithms is the selection process. The NSGA-III uses reference points to select solutions to the next generation.

input_data folder:

This folder contains the definition of all problems. The format is also extended form Kanpur Genetic Algorithms Laboratory.

For example: car side impact problem definition.

2                -->:  Solver Type (2 --> A^2 NSGA-III, 1 --> A NSGA-III, 0 --> NSGA-III)

156              -->:  Population Size (a multiple of 4, because it is necessary for tournament selection routine)
500              -->:  Number of generations
3                -->:  Number of objectives
10               -->:  Number of inequality constrains
0                -->:  Number of equality constrains
7                -->:  Number of real variables
0.5 1.5          -->:  Low limit to high limit of real variable 1
0.45 1.35        -->:  Low limit to high limit of real variable 2
0.5 1.5          -->:  Low limit to high limit of real variable 3
0.5 1.5          -->:  Low limit to high limit of real variable 4
0.875 2.625      -->:  Low limit to high limit of real variable 5
0.4 1.2          -->:  Low limit to high limit of real variable 6
0.4 1.2          -->:  Low limit to high limit of real variable 7
0.1              -->:  Crossover probability of real variable (0.6-1.0)
1                -->:  Mutation probablity of real variables (1/nreal)
20               -->:  Distribution index for crossover (5-20)
30               -->:  Distribution index for mutation (5-50)
0                -->:  Number of binary variables
1                -->:  Do you want to use gnuplot to display the results realtime (0 - NO) (1 - yes)
3                -->:  2D display, 3D display or Parallel Coordinates?, enter 2 for 2D, 3 for 3D and P for Parallel coordinates
1                -->:  Objective for X display
2                -->:  Objective for y display
3                -->:  Objective for Z display
30               -->:  Enter the first angle (an integer in the range 0-180)
60               -->:  Enter the first angle (an integer in the range 0-360)

real_front folder:
