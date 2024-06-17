#ifndef POWERMETHOD_HPP
#define POWERMETHOD_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <ctime>


using namespace std;

class PowerMethod
{
private:
    float** M;
    float** U;
    float* x;
    int dim;
public:
/// Constructor
    PowerMethod();
    PowerMethod(string path);
    PowerMethod(int dim, float** M);
/// Destructor
    ~PowerMethod();
/// Display
    void  display(int precision = 3);
    void solution(int precision = 3);
/// Fill data
    void fill(string path);
    void fill(int dim, float** M);
/// Operation
    void dot(float* vec, float* buff);
    float product(float* v1, float* v2);
/// Swap
    void swap(float** v1, float** v2);
/// Normalization
    float norm(float* vec);
    void  normalize(float* vec);
/// Find Eigenvalues
    void compute(int k);
    void deflate();
};

#endif