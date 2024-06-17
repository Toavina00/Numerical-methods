#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>


using namespace std;

template <typename T>
class Solver
{
private:
    T** M;
    T *b, *x, *y;
    int dim;
    bool resolved;
    bool factored;
public:
/// Constructor
    Solver();
    Solver(string path);
    Solver(int dim, T** M, T* b);
/// Destructor
    ~Solver();
/// Display
    void display(int precision = 3);
    void solution(int precision = 3);
/// Fill data
    void fill(string path);
    void fill(int dim, T** M, T* b);
/// Swap rows
    void swap(int r0, int r1);
/// Solve
    void factorization();
    void solve();
};

#endif