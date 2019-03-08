# NSGA-III- A Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems (Non-Official)
This work provides a third evolution Nondominated Sorting Genetic Algorithm (NSGA-III) implementation, extending the official NSGA-II algorithm stored in the repository of the Kanpur Genetic Algorithms Laboratory (KanGAL) and modifying the selection operator. The algorithm resolves binary and real, constrained and non-constrained Evolutionary Multi-Objective and Many-objective Optimization problems. 
Efficiently Adaptive (A^2-NSGA-III) and adaptive (A-NSGA-III) third evolution Nondominated Sorting Genetic Algorithms are extended as well. They solve the non-uniform random Pareto-optimal distribution inconvenient of some standardized tests and probably most of real problems, for which some useless reference points distributed on the M dimensional hyperplane, are not associated to any population member. Ideally, we desire the niche count of all reference points is equal to 1, to obtain a uniform Pareto-optimal Front. Also, the adaptive Reference points generation for A-NSGA-III and A^2-NSGA-III algorithms is expanded. Our NSGA-III algorithm outperforms most of the author results for standardized DTLZ problems in terms of the Inverted Generational Distance measure. In addition, real problems such as the Car-Side Impact and Water problems find a better Pareto-optimal distribution, not only visually but also in terms of the hypervolume measure, when adaptive and efficiently adaptive reference points are employed.

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
