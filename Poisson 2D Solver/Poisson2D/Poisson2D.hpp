#ifndef POISSON2D_HPP
#define POISSON2D_HPP

#include <functional>
#include <cmath>
#include <iomanip>
#include "Tools.hpp"


/*******************************/
/****** Poisson 2D Solver ******/
/*******************************/

/*
    Solve the following problem
        -u" = f
        on a domain [a,L]x[b,W] with boundary conditions:
            u(a,y) = g(a,y)
            u(L,y) = g(L,y)
            u(x,b) = g(x,b)
            u(x,W) = g(x,W)
*/

using namespace std;

class Poisson2D
{
private:
    double x, y;                     // Discretisation grid origin
    double x_len, y_len;             // Discretisation grid size
    double dx, dy;                   // Discretisation steps
    int N, M;                        // Discretisation grid subdivision number
    
public:
    Poisson2D();
    ~Poisson2D();

    void setDomain(double x, double x_len, double y, double y_len, int M, int N);               // Define the rectangular domain
    function<double(int, int)> discretize(function<double(double, double)> func);               // Discretize  the given function on a regular grid
    
    double** solve(function<double(double, double)> f, function<double(double, double)> g);     // Solve the problem given f and the boundary condition
};


#endif