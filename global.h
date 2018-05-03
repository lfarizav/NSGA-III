/* This file contains the variable and function declarations */

# ifndef _GLOBAL_H_
# define _GLOBAL_H_

# define INF 1.0e14
# define EPS 1.0e-14
# define E  2.71828182845905
# define PI 3.14159265358979
# define GNUPLOT_COMMAND "gnuplot -persist"

typedef struct
{
    int size_refpoints;
    int size_refpoints_inside;
    int size_useless_refpoints_number;
    int size_usefull_refpoints_number;
    int points_per_dim_number;
    double step;
}
individual_refpoints;

typedef struct
{
    individual_refpoints *ind;
}
population_refpoints;

typedef struct
{
    int rank;
    double constr_violation;
    double equality_constr_violation;
    double *xreal;
    double *xbin;
    int **gene;
    double *obj;
    double *obj_minus_zmin;
    double *obj_normalized;
    double *obj_feasible;
    double *obj_infeasible;
    double *constr;
    double *equality_constr;
    int associatedref;
    double distancetoassociatedref;
    double w;
    int is_feasible;
}
individual;

typedef struct
{
    individual *ind;
}
population;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
}
list;

extern int nreal;
extern int nbin;
extern int nobj;
extern int ncon;
extern int neqcon;
extern int popsize;
extern int adaptive_nsga;
extern double pcross_real;
extern double pcross_bin;
extern double pmut_real;
extern double pmut_bin;
extern double eta_c;
extern double eta_m;
extern double *a_last_gen;
extern int ngen;
extern int nbinmut;
extern int nrealmut;
extern int nbincross;
extern int nrealcross;
extern int *nbits;
extern double *min_realvar;
extern double *max_realvar;
extern double *min_binvar;
extern double *max_binvar;
extern int bitlength;
extern int choice;
extern int obj1;
extern int obj2;
extern int obj3;
extern int obj4;
extern int obj5;
extern int angle1;
extern int angle2;
extern double *scale_obj_min;
extern double *scale_obj_max;
extern int *scale_obj_min_ref;
extern int *scale_obj_max_ref;
extern int dtlz;
extern double *a;
extern double *smin;
extern int *index_s;
extern double **zmax;
extern int numberpointperdim;
extern int numberpointperdim_inside;
extern long factorial;
extern long factorial_inside;
extern long factorial_adaptive;
extern int *rho;
extern int *rho_St;
extern int *rho_Fl;
extern int *last_rho_St;
extern int *last_generation_associated_rho_St;
extern int min_rho_St;
extern int min_rho_Fl;
extern int min_rho;
extern int min_rho_Fl_count_index;
extern int min_rho_count_index;
extern int *associatedfromlastfront_St;
extern int *associatedfromlastfront_Fl;
extern int *membertoadd;
extern int *ref_points_min_rho;
extern int *ref_points_min_rho_Fl;
extern FILE gp_test;
extern double *w_scalarizing_vector;
extern int first_front;
extern int *fronts;
extern double **ref_points;
extern double **adaptive_ref_points_settled;
extern int *adaptive_ref_points_settled_number;
extern int elegible_adaptive_ref_points_to_be_fixed_number;
extern double **DTLZ;
extern double **temp_refpoints;
extern int *temp_refpoints_pointer;
extern double **minimum_amount_ref_points;
extern double **ref_points_normalized;
extern double *min_ref_points;
extern double *max_ref_points;
extern int *dist_lf;
extern int *index;
extern double **igb_real_front;
extern double **igb_algorithm;
extern double **igb_algorithm_normalized;
extern double **igb_real_front_normalized;
extern double *maximum_value;
extern double *minimum_value;
extern int number_is_feasible;
extern int number_is_infeasible;
extern int *infeasible_population_sorted_list_index;
extern int *feasible_population_sorted_list_index;
extern int *useless_refpoint_index;
extern int *usefull_refpoint_index;
extern int *sort_all_refpoint_index;
extern int *sort_all_adaptive_refpoint_index;
extern int useless_refpoint_number;
extern int usefull_refpoint_number;
extern double **adaptive_refpoints;
extern int adaptive_refpoint_number;
extern int last_gen_adaptive_refpoints_number;
extern double *num_div_den;
extern int IGDfrontsize;
extern int counter;
extern int adaptive_ref_points_inserted;
extern int adaptive_ref_points_inserted_per_generation;
void allocate_memory_pop (population *pop, int size);
void allocate_memory_ind (individual *ind);
void deallocate_memory_pop (population *pop, int size);
void deallocate_memory_ind (individual *ind);

double maximum (double a, double b);
double minimum (double a, double b);

void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2);
void realcross (individual *parent1, individual *parent2, individual *child1, individual *child2);
void bincross (individual *parent1, individual *parent2, individual *child1, individual *child2);

void decode_pop (population *pop);
void decode_ind (individual *ind);

void onthefly_display (population *pop, FILE *gp, int ii, int normalization);
void onthefly_display_real_front (population *pop,FILE *gp);
void onthefly_display_minus_zmin (population *pop, FILE *gp);
void onthefly_display_normalized (population *pop, FILE *gp, int archieve_and_front_sizes);
void onthefly_display_a (population *pop, FILE *gp);
void onthefly_display_one (population *pop, FILE *gp);
void onthefly_display_inside (population *pop, FILE *gp);
void onthefly_display_zmax (population *pop, FILE *gp);
void onthefly_display_DTLZ1 (FILE *gp_dtlz);
void onthefly_display_refpoints (population *pop, FILE *gp);
void onthefly_display_parallel_coordinates (population *pop, FILE *gp_pc, int ii);
void display_pop (population *pop);
void display_pop_ind_obj (individual *ind, int popsizeline);
void display_pop_ind_obj_minus_zmin (individual *ind, int popsizeline);
void display_pop_ind_obj_normalized (individual *ind, int popsizeline);
void display_pop_ind_xreal (individual *ind, int popsizeline);
void display_refpoints ();
void display_DTLZ1 ();
void display_refpoints_normalized ();
void display_fronts ();
void display_d (int pop_size);

int check_dominance (individual *a, individual *b);

void evaluate_pop (population *pop);
void evaluate_ind (individual *ind);

void fill_nondominated_sort (population *selection_pop, population *mixed_pop, population *new_pop, int generation);
int bubble_sorting_infeasible_population_index(population *poputation_sorted);
void feasible_population_index(population *poputation_sorted);
/*new selection for NSGA-III*/
void associated_reference_points_fill (population *selection_pop, population *mixed_pop, population *new_pop, int front_size,int achieve_size, list *cur, int generation);

void reset_pop_obj (population *pop,int size);
void reset_ind_obj (individual *ind);
void initialize_pop (population *pop);
void initialize_ind (individual *ind);

void insert (list *node, int x);
list* del (list *node);

void merge(population *pop1, population *pop2, population *pop3);
void copy_ind (individual *ind1, individual *ind2);

void mutation_pop (population *pop);
void mutation_ind (individual *ind);
void bin_mutate_ind (individual *ind);
void real_mutate_ind (individual *ind);

void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, double *equality_constr, int normalized);

void report_pop (population *pop, FILE *fpt);
void report_feasible (population *pop, FILE *fpt);
void report_ind (individual *ind, FILE *fpt);

void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size);
void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right);
void quicksort_dist(population *pop, int *dist, int front_size);
void q_sort_dist(population *pop, int *dist, int left, int right);

void selection (population *old_pop, population *new_pop);
individual* tournament (individual *ind1, individual *ind2);

/*selection.c- NSGA-III*/
int variables_initialization(population *selection_pop, population *mixed_pop, population *new_pop, int front_size, int archieve_size, list *elite);
void find_min_from_functions(individual *ind, int i, int population_type);
void find_max_from_functions(individual *ind, int i, int population_type);
void obj_minus_zmin (individual *ind);
void find_a ();
void cofactor(double num[25][25], int f);
void transpose(double num[25][25], double fac[25][25], int r);
double determinant(double zmax_matrix[25][25], int nobj);
void normalized_objective_function (individual *ind);
void normalized_objective_function_simple (individual *ind);
int niching (population *selection_pop,population *new_pop, int front_size, int archieve_size, int start, int end);
int associated_from_last_front(individual *normalizedind, int l, int index, int associatedfromlastfront_index);
void associate(individual *normalizedind,individual *new_ind,int l, int archieve_size, int start, int end);
void get_scalarizing_vector(int j);
double achievement_scalarization_function (individual *ind_minus_zmin, int i);
void find_extreme_points(population *selection_pop_minus_zmin, int archieve_size);
void find_max_objectives(individual *ind, int i);
void construct_hyperplane ();
int is_zmax_duplicated ();
double perpendicular_distance(individual *normalizedind, int l);

void invert_real_front(int dtlz);
void get_fronts_from_file (int dtlz);
double IGD (population *pop);
void get_maximum_value(int front);
void get_minimum_value(int front);
void get_normalized_front (int front);

int generate_DTLZ1 (int nobj_for, double step);
int recursive_for_DTLZ1 (int nobj_for, double step, int count, int i);
int generate_ref_points (int nobj_for, double step);
int recursive_for (int nobj_for, double step, int count, int i);
int generate_ref_points_inside (int nobj_for, double step);
int recursive_for_inside (int nobj_for, double step, int count, int i);
int generate_adaptive_ref_points_inside (int nobj_for, double step);
int recursive_for_adaptive (int nobj_for, double step, int count, int i);
int create_adaptive_refpoints();
long fact(int x);
void refpoints_normalized();
void min_refpoints();
void max_refpoints();
void square_refpoints();
int get_supplied_reference_points_from_file (char *filename);
void sort_all_refpoints_by_rho_index();
void add_adaptive_refpoints_to_ref_points();
int delete_adaptive_refpoints(int archieve_size, int front_size, population *selection_pop, population *new_pop, int generation);
void find_useless_usefull_refpoint_index();
void check_adaptive_refpoints_inclusion_number(int generation);
void store_useless_refpoints (int adaptive_ref_point_number, int useless_ref_point_number);
void load_useless_refpoints(int adaptive_ref_point_number, int useless_ref_point_number);
# endif
