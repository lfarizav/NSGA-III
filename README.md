# NSGA-III: A Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems (Non-Official)
This work provides a third evolution Nondominated Sorting Genetic Algorithm (NSGA-III) implementation, extending the official NSGA-II algorithm stored in the repository of the Kanpur Genetic Algorithms Laboratory (KanGAL) and modifying the selection operator. The algorithm resolves binary and real, constrained and non-constrained Evolutionary Multi-Objective and Many-objective Optimization problems. 
Efficiently Adaptive (A^2-NSGA-III) and adaptive (A-NSGA-III) third evolution Nondominated Sorting Genetic Algorithms are extended as well. They solve the non-uniform random Pareto-optimal distribution inconvenient of some standardized tests and probably most of real problems, for which some useless reference points distributed on the M dimensional hyperplane, are not associated to any population member. Ideally, we desire the niche count of all reference points to be equal to 1, obtaining a uniform Pareto-optimal Front. Also, the adaptive Reference points generation for A-NSGA-III and A^2-NSGA-III algorithms is expanded. Our NSGA-III algorithm outperforms most of results for standardized DTLZ problems in terms of the Inverted Generational Distance measure. In addition, real problems such as the Car-Side Impact and Water problems find a better Pareto-optimal distribution, not only visually but also in terms of the hypervolume measure, when adaptive and efficiently adaptive reference points are employed.
## Introducction
Genetic algorithms (GAs) are population-based stochastic search and optimization of evolutionary computation methods that resolve problems in areas such as Engineering, Sciences, and Commerce \cite{introgenalg1999}. GAs are preferred from Classical optimization methods (Direct and Gradient-based), due to both, the versatility to resolve different complex problems and the facility to find multiple solutions without a priori problem information. They are inspired by genetic fundamental laws such as the natural selection, first introduced by John Holland and his colleagues at Michigan University in the 60's and 70's. Further, Goldberg developed and popularized Holland's methods in 1989 \cite{compintelligence}. 
The basic framework of GAs considers an initial population of pseudo-random solutions. Then, they pass through genetic functions such as selection, crossover, and mutation to recombine and perturb solutions. Finally, solutions are evaluated using a \textit{fitness function} in the hope of creating better solutions. The fittest solutions survive and evolve to the next generation, but the remaining are neglected. This process is repeated until a predefined termination criterion is met.

Many real-life problems have more than one conflicting objective functions to be optimized, so either a multi-objective optimization problem (MOOP) or a many-objective optimization problem (MaOP) formulation is available.  Solutions must promote solutions diversity, use specialized fitness functions for evaluations and satisfy the objective functions within an acceptable degree without being dominated by other solutions \cite{GAtutorial}. 

Multi-Objective Genetic Algorithms, use a modified single objective GA to find multiple nondominated solutions. Multi-Objective Genetic Algorithms have the property of searching simultaneously different regions of non-convex, discontinuous, and multimodal solution spaces.  Some Multi-Objective Genetic Algorithms are: Vector Evaluated Genetic Algorithm (VEGA) \cite{vega}, Multi-Objective Genetic Algorithms (MOGA), Niched Pareto Genetic Algorithm (NPGA), Weight-Based Genetic Algorithm (WBGA), Random Weighted Genetic Algorithm (RWGA), Nondominated Sorting Genetic Algorithm (NSGA) \cite{nsga}, Strength Pareto Evolutionary Algorithm (SPEA), Pareto Archived Evolution Strategy (PAES), Pareto Enveloped-Based Selection Algorithm (PESA), Multi-objective Evolutionary Algorithm (MEA), Micro GA, Rank-Density Based Genetic Algorithm (RDGA), Dynamic Multi-Objective Evolutionary Algorithm (DMOEA) \cite{GAtutorial}. Further improvements are the Cellular MOGA, SPEA2 \cite{spea2}, PESA-II \cite{pesa2}, NSGA-II, \cite{nsgaii} and NSGA-III \cite{nsgaiiipart1}\cite{nsgaiiipart2}. Some available on-line implementations from scratch without using Deb's official NSGA-II code, are the Yarpiz (Matlab) \cite{yarpiz}, jMetal (Java) \cite{jmetal}, nsga3cpp (C++) \cite{nsga3cpp} and nsga3 (Python) \cite{nsga3python}.
Typically, the performance evaluation for evolutionary MaOPs is measured using metrics like the Generational Distance (GD), the Inverted Generational Distance (IGD), or the hypervolume to provide a comprehensive comparison among algorithms.
## Motivation
A wide range of real-world constrained and unconstrained problems involve more than 3 dimensions. MOOPs well-established algorithms are not suited when dealing with this kind of problems, basically because of the dimensionality curse, therefore it behaves similarly to a random walk. Nevertheless, there is a MaOP formulation that tackles such high-dimensional optimization issues.

The NSGA-III is a MaOP algorithm that uses well-distributed reference points, extending the NSGA-II selection method. It tries to overcome the increase of nondominated solutions with the number of objectives functions, ineffective mutation and crossover operations, the diversity measure estimation problem and the high-dimensional visualization \cite{Springer:RAEMOptimization}. It has the ability to find a well-converged and well-diversified set of solutions in high-dimensional scenarios. Additionally, it reduces the computational complexity, and increases the efficiency up to 15 objective functions.
However, the official code is proprietary \cite{nsga3proprietary} and it is not shared with the research community. It is hard and inconvenient for people who want to either apply or extend MaOP algorithms to real problems.

I decided to use the official NSGA-II KanGAL code to build the NSGA-III, A-NSGA-III and A^2-NSGA-III versions, because I think it is the smoothest way. Other implementations shared on the Internet start from scratch which increase the complexity level. Only the reference \cite{nsga3cpp} shows validation results for DTLZ problems. Furthermore, I determine to complement the explanation about adaptive reference points generation.
The intention of this paper is to build the most similar NSGA-III code compared with the official and share the code to the research community. I would like to apply it to resolve Cloud Radio Access Network (C-RAN) planning problems to reduce the power consumption in cellular networks. 
I will expect this paper contribute to many researchers to solve many questions with regard to the implementation of the original NSGA-III. 
## Instructions
### input_data folder:

This folder contains the definition of all problems. The Kanpur Genetic Algorithms Laboratory format form is slighly modified.

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
