
// Generate Pareto front for DTLZ1-4 functions following NSGA-III's paper.
//
// An Evolutionary Many-Objective Optimization Algorithm Using Reference-point Based Non-dominated Sorting Approach, Part I:
// Solving Problems with Box Constraints
//
// http://dx.doi.org/10.1109/TEVC.2013.2281535

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

typedef vector<double> TObjVec;
typedef vector<TObjVec> TFront;

void GeneratePF_OneLayer(ostream &ofile, const string &problem_name, int M, int p);
void GeneratePF_TwoLayers(ostream &os, const string &problem_name, int M, int outside_p, int inside_p);
string IntToStr(int i);

int main()
{
    const char *problem_name[]= {"DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4"};

    for (int i=0; i<4; i+=1) // problem
    {
        const int M[5] = {3, 5, 8, 10, 15};
        for (int j=0; j<5; j+=1) // objectives
        {
            ofstream ofile( string(problem_name[i]) + "(" + IntToStr(M[j]) + ")-PF.txt");

            if (M[j] <= 5) // #objectives <= 5
            {
                int p[2] = {12, 6}; // Check Section V, Table I in the original paper
                GeneratePF_OneLayer(ofile, problem_name[i], M[j], p[j]);
            }
            else
            {
                int p[3][2] = {{3, 2}, {3, 2}, {2, 1}}; // Check Section V, Table I in the original paper
                GeneratePF_TwoLayers(ofile, problem_name[i], M[j], p[j-2][0], p[j-2][1]);
            }
        }
    }

    return 0;
}

// ----------------------------------------------------------------------
void generate_recursive(TFront *pf, TObjVec *pt, size_t num_objs,
                        size_t left, size_t total, size_t element)
{
    if (element == num_objs-1)
    {
        (*pt)[element] = left;
        pf->push_back(*pt);
    }
    else
    {
        for (size_t i=0; i<=left; i+=1)
        {
            (*pt)[element] = i;
            generate_recursive(pf, pt, num_objs, left-i, total, element+1);
        }
    }
}
// ----------------------------------------------------------------------
void GenerateWeight(TFront *pf, size_t M, size_t p)
{
    TObjVec pt(M);

    generate_recursive(pf, &pt, M, p, p, 0);
}
// ----------------------------------------------------------------------
void GeneratePF_OneLayer(ostream &os, const string &problem_name, int M, int p)
{
    TFront PF;

    int num_objectives = M, num_divisions = p;
    GenerateWeight(&PF, num_objectives, num_divisions);

    if (problem_name == "DTLZ1")
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << (0.5*PF[i][j])/num_divisions << ' ';
            }
            os << endl;
        }
    }
    else // DTLZ2-4
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            double sum = 0;

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                sum += PF[i][j]*PF[i][j];
            }

            double k = sqrt(1.0/sum);

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << k*PF[i][j] << ' ';
            }
            os << endl;
        }
    } // else

}// GeneratePF_OneLayer()
// ----------------------------------------------------------------------
void GeneratePF_TwoLayers(ostream &os, const string &problem_name, int M, int outside_p, int inside_p)
{

    GeneratePF_OneLayer(os, problem_name, M, outside_p);


    TFront PF;

    int num_objectives = M, num_divisions = inside_p;
    GenerateWeight(&PF, num_objectives, num_divisions);

    for (size_t i=0; i<PF.size(); i+=1)
    {
        for (size_t j=0; j<PF[i].size(); j+=1)
        {
            PF[i][j] = (static_cast<double>(num_divisions)/M+PF[i][j])/2; // (k=num_divisions/M, k, k, ..., k) is the center point
        }
    }

    if (problem_name == "DTLZ1")
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << (0.5*PF[i][j])/num_divisions << ' ';
            }
            os << endl;
        }
    }
    else // DTLZ2-4
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            double sum = 0;

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                sum += PF[i][j]*PF[i][j];
            }

            double k = sqrt(1.0/sum);

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << k*PF[i][j] << ' ';
            }
            os << endl;
        }
    } // else

}// GeneratePF_TwoLayers()
// ----------------------------------------------------------------------

string IntToStr(int i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}
