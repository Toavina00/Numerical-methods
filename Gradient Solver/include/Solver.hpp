#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <random>

using namespace std;

template <typename T>
class Solver
{
private:
    T** M;
    T *b, *x, *y, *r;
    int dim;
    bool resolved;
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
/// Operations
    void add(T* a, T* b, T* buff);
    void sub(T* a, T* b, T* buff);
    void scalar(T sc, T* a, T* buff);
    void dot(T* a, T* buff);
    T product(T* a, T* b);
    T norm(T* a);
/// Solve
    void solve();
};

#endif