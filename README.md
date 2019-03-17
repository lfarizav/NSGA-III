## Copyright
> NSGA-III and NSGA-II Copyrights belong to Luis Felipe Ariza Vesga and the Kanpur Genetic Algorithms Laboratory, respectively. You are free to use this algorithm (https://github.com/lfarizav/NSGA-III) for **research purposes**. All publications which use this code should **acknowledge the author**. Luis Felipe Ariza Vesga. A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems. March, 2019.
```
@online{NonofficialNSGAIII,
  title={A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems},
  author={Luis Felipe Ariza Vesga},
  url = {https://github.com/lfarizav/NSGA-III}
  month = March,
  year={2019},
}
```
# NSGA-III: A Fast Nondominated Sorting Genetic Algorithm Extension to Solve Evolutionary Many-Objective Problems (Non-Official)
This work provides a third evolution fast Nondominated Sorting Genetic Algorithm (NSGA-III) implementation in C, extending the official NSGA-II algorithm stored in the repository of the Kanpur Genetic Algorithms Laboratory (KanGAL) (https://www.iitk.ac.in/kangal/codes.shtml). The algorithm resolves binary and real, constrained and non-constrained Evolutionary Multi-Objective and Many-objective Optimization problems. 
Efficiently Adaptive (A^2-NSGA-III) and adaptive (A-NSGA-III) third evolution Nondominated Sorting Genetic Algorithms are extended as well. They solve the non-uniform random Pareto-optimal distribution inconvenient of some standardized tests and probably most of real problems, for which some useless reference points distributed on the M dimensional hyperplane, are not associated to any population member. Ideally, we desire the niche count of all reference points to be equal to 1, obtaining a uniform Pareto-optimal Front. My NSGA-III algorithm outperforms most of results for standardized DTLZ problems in terms of the Inverted Generational Distance measure. In addition, real problems such as the Car-Side Impact and Water problems find a better Pareto-optimal distribution, not only visually but also in terms of the hypervolume measure, when adaptive and efficiently adaptive reference points are employed.
## Introducction
Genetic algorithms (GAs) are population-based stochastic search and optimization of evolutionary computation methods that resolve problems in areas such as Engineering, Sciences, and Commerce \cite{introgenalg1999}. GAs are preferred from Classical optimization methods (Direct and Gradient-based), due to both, the versatility to resolve different complex problems and the facility to find multiple solutions without a priori problem information. They are inspired by genetic fundamental laws such as the natural selection, first introduced by John Holland and his colleagues at Michigan University in the 60's and 70's. Further, Goldberg developed and popularized Holland's methods in 1989 \cite{compintelligence}. 
The basic framework of GAs considers an initial population of pseudo-random solutions. Then, they pass through genetic functions such as selection, crossover, and mutation to recombine and perturb solutions. Finally, solutions are evaluated using a fitness function in the hope of creating better solutions. The fittest solutions survive and evolve to the next generation, but the remaining are neglected. This process is repeated until a predefined termination criterion is met.

Many real-life problems have more than one conflicting objective functions to be optimized, so either a multi-objective optimization problem (MOOP) or a many-objective optimization problem (MaOP) formulation is available.  Solutions must promote solutions diversity, use specialized fitness functions for evaluations and satisfy the objective functions within an acceptable degree without being dominated by other solutions \cite{GAtutorial}. 

Multi-Objective Genetic Algorithms, use a modified single objective GA to find multiple nondominated solutions. Multi-Objective Genetic Algorithms have the property of searching simultaneously different regions of non-convex, discontinuous, and multimodal solution spaces.  Some Multi-Objective Genetic Algorithms are: Vector Evaluated Genetic Algorithm (VEGA) \cite{vega}, Multi-Objective Genetic Algorithms (MOGA), Niched Pareto Genetic Algorithm (NPGA), Weight-Based Genetic Algorithm (WBGA), Random Weighted Genetic Algorithm (RWGA), Nondominated Sorting Genetic Algorithm (NSGA) \cite{nsga}, Strength Pareto Evolutionary Algorithm (SPEA), Pareto Archived Evolution Strategy (PAES), Pareto Enveloped-Based Selection Algorithm (PESA), Multi-objective Evolutionary Algorithm (MEA), Micro GA, Rank-Density Based Genetic Algorithm (RDGA), Dynamic Multi-Objective Evolutionary Algorithm (DMOEA) \cite{GAtutorial}. Further improvements are the Cellular MOGA, SPEA2 \cite{spea2}, PESA-II \cite{pesa2}, NSGA-II, \cite{nsgaii} and NSGA-III \cite{nsgaiiipart1}\cite{nsgaiiipart2}. Some available on-line NSGA-III implementations from scratch without using Deb's official NSGA-II code, are the Yarpiz (Matlab) \cite{yarpiz}, jMetal (Java) \cite{jmetal}, nsga3cpp (C++) \cite{nsga3cpp}, nsga3 (Python) \cite{nsga3python} and Platemo [24].
Typically, the performance evaluation for evolutionary MaOPs is measured using metrics like the Generational Distance (GD), the Inverted Generational Distance (IGD), or the hypervolume to provide a comprehensive comparison among algorithms.
## Motivation
A wide range of real-world constrained and unconstrained problems involve more than 3 dimensions. MOOPs well-established algorithms are not suited when dealing with this kind of problems, basically because of the dimensionality curse, therefore it behaves similarly to a random walk. Nevertheless, there is a MaOP formulation that tackles such high-dimensional optimization issues.

The NSGA-III is a MaOP algorithm that uses well-distributed reference points, extending the NSGA-II selection method. It tries to overcome the increase of nondominated solutions with the number of objectives functions, ineffective mutation and crossover operations, the diversity measure estimation problem and the high-dimensional visualization \cite{Springer:RAEMOptimization}. It has the ability to find a well-converged and well-diversified set of solutions in high-dimensional scenarios. Additionally, it reduces the computational complexity, and increases the efficiency up to 15 objective functions.

However, the official code is proprietary and it is not shared with the research community. It is hard and inconvenient for people who want to either apply or extend MaOP algorithms to real problems.

I decided to use the official NSGA-II KanGAL code to build the NSGA-III, A-NSGA-III and A^2-NSGA-III versions, because I think it is the smoothest way and it is implemented in C. Other shared implementations in Matlab are not fast, so the C solution is preferred.

The intention of this paper is to build the most similar NSGA-III code compared with the official and share the code to the research community. I would like to apply it to resolve Cloud Radio Access Network (C-RAN) planning problems to reduce the power consumption in cellular networks. 

I will expect this paper contributes to many researchers to solve many questions with regard to the implementation of the original NSGA-III. 

## Instructions
### input_data:

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

### real_front:
This folder has as much as possible dtlz fronts.
### Hou to use:
> make nsga2r

> ./nsga2r seed <input_data/dtlz1.in ----> for the dtlz1 problem.

## References
[1] K. Deb, “An introduction to genetic algorithms,” Sadhana, vol. 24, no. 4, pp. 293–315, Aug 1999. [Online]. Available: https://doi.org/10. 1007/BF02823145

[2] J. Kacprzyk and W. Pedrycz, Springer Handbook of Computational Intelligence. Springer Publishing Company, Incorporated, 2015.

[3] A. Konak, D. W. Coit, and A. E. Smith, “Multi-objective optimization using genetic algorithms: A tutorial,” Reliability Engineering & System Safety, vol. 91, no. 9, pp. 992 – 1007, 2006, special Issue - Genetic Algorithms and Reliability. [Online]. Available: http://www.sciencedirect.com/science/article/pii/S0951832005002012

[4] J. D. Schaffer, “Multiple objective optimization with vector evaluated genetic algorithms,” in Proceedings of the 1st International Conference on Genetic Algorithms. Hillsdale, NJ, USA: L. Erlbaum Associates Inc., 1985, pp. 93–100. [Online]. Available: http://dl.acm.org/citation. cfm?id=645511.657079

[5] N. Srinivas and K. Deb, “Muiltiobjective optimization using nondominated sorting in genetic algorithms,” Evolutionary Computation, vol. 2, no. 3, pp. 221–248, Sept 1994.

[6] M. Kim, T. Hiroyasu, M. Miki, and S. Watanabe, “Spea2+: Improving the performance of the strength pareto evolutionary algorithm 2,” in Parallel Problem Solving from Nature - PPSN VIII, X. Yao, E. K. Burke, J. A. Lozano, J. Smith, J. J. Merelo-Guervos, J. A. Bullinaria, J. E. ´ Rowe, P. Tino, A. Kab ˇ an, and H.-P. Schwefel, Eds. Berlin, Heidelberg: ´ Springer Berlin Heidelberg, 2004, pp. 742–751.

[7] D. W. Corne, N. R. Jerram, J. D. Knowles, and M. J. Oates, “Pesa-ii: Region-based selection in evolutionary multiobjective optimization,” in Proceedings of the 3rd Annual Conference on Genetic and Evolutionary Computation, ser. GECCO’01. San Francisco, CA, USA: Morgan Kaufmann Publishers Inc., 2001, pp. 283–290. [Online]. Available: http://dl.acm.org/citation.cfm?id=2955239.2955289

[8] N. Srinivas and K. Deb, “Muiltiobjective optimization using nondominated sorting in genetic algorithms,” Evolutionary Computation, vol. 2, no. 3, pp. 221–248, Sept 1994.

[9] K. Deb and H. Jain, “An evolutionary many-objective optimization algorithm using reference-point-based nondominated sorting approach, part i: Solving problems with box constraints,” IEEE Transactions on Evolutionary Computation, vol. 18, no. 4, pp. 577–601, Aug 2014.

[10] H. Jain and K. Deb, “An evolutionary many-objective optimization algorithm using reference-point based nondominated sorting approach, part ii: Handling constraints and extending to an adaptive approach,” IEEE Transactions on Evolutionary Computation, vol. 18, no. 4, pp. 602–622, Aug 2014.

[11] Yarpiz. (2018) Nsga-iii: Non-dominated sorting genetic algorithm, the third version. [Online]. Available: http://yarpiz.com/456/ypea126-nsga3

[12] jMetal. (2018) Metaheuristic algorithms in java. [Online]. Available: http://jmetal.sourceforge.net/

[13] T.-C. Chiang and Collaborators. (2014) nsga3cpp: A c++ implementation of nsga-iii. [Online]. Available: http://web.ntnu.edu.tw/∼tcchiang/publications/nsga3cpp/nsga3cpp.htm

[14] L. Marti. (2016) nsga3cpp: A c++ implementation of nsga-iii. [Online]. Available: https://github.com/lmarti/nsgaiii

[15] L. B. S. Slim Bechikh, Maha Elarbi, A Survey in Recent Advances in Evolutionary Multi-objective Optimization, ser. Many-objective Optimization Using Evolutionary Algorithms:. New York, NY: Springer, 2017, vol. 20. [Online]. Available: http://dx.doi.org/8843/10. 1007/978-3-319-42978-6

[16] R. C. T. (Siemens). (2018) Heeds smashes barriers on multi-objective design studies. [Online]. Available: https://www.redcedartech.com/ newsletters/HEEDS News-Mar15.htm

[17] Y. Yusoff, M. S. Ngadiman, and A. M. Zain, “Overview of nsga-ii for optimizing machining process parameters,” Procedia Engineering, vol. 15, pp. 3978 – 3983, 2011, cEIS 2011. [Online]. Available: http://www.sciencedirect.com/science/article/pii/S1877705811022466

[18] I. Das and J. E. Dennis, “Normal-boundary intersection: A new method for generating the pareto surface in nonlinear multicriteria optimization problems,” SIAM J. on Optimization, vol. 8, no. 3, pp. 631–657, Mar. 1998. [Online]. Available: http://dx.doi.org/10.1137/
S1052623496307510

[19] H. Jain and K. Deb, “An improved adaptive approach for elitist nondominated sorting genetic algorithm for many-objective optimization,” in Evolutionary Multi-Criterion Optimization, R. C. Purshouse, P. J. Fleming, C. M. Fonseca, S. Greco, and J. Shaw, Eds. Berlin, Heidelberg: Springer Berlin Heidelberg, 2013, pp. 307–321.

[20] S. Jiang and S. Yang, “A strength pareto evolutionary algorithm based on reference direction for multiobjective and many-objective optimization,” IEEE Transactions on Evolutionary Computation, vol. 21, no. 3, pp.
329–346, June 2017.

[21] K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, Scalable Test Problems for Evolutionary Multiobjective Optimization. London: Springer London, 2005, pp. 105–145. [Online]. Available: https: //doi.org/10.1007/1-84628-137-7 6

[22] S. Z. P. S. W. L. Q.Zhang, A. Zhou and S. Tiwari, “Multiobjective optimization test instances for the cec-2009 special session and competition,” Nanyang Technol. Univ., Singapore, Tech, 2008. [Online]. Available: http://www.nyu.edu.sg/home/epnsugan

[23] C. M. Fonseca, L. Paquete, and M. Lopez-Ib ´ a´nez, “An improved dimension-sweep algorithm for the hypervolume indicator,” in Proceedings of the 2006 Congress on EvolutionaryComputation (CEC 2006). Piscataway, NJ: IEEE Press, July 2006, pp. 1157–1163.

[24] Y. Tian, R. Cheng, X. Zhang, and Y. Jin, “Platemo: A matlab platform for evolutionary multi-objective optimization [educational forum],” IEEE Computational Intelligence Magazine, vol. 12, no. 4, pp. 73–87, Nov 2017.

