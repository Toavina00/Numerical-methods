#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

using namespace std;


/****** Display functions ******/
void display(double* x, int N);
void display(double** x, int N, int M);
void display(complex<double>* x, int N);

/****** Array operations ******/
template <typename T> void       factor(T* X, double x, int N);

/****** FFT and IFFT ******/
template <typename T> void   fourier_ft(complex<double>* X, T* x, int N);
template <typename T> void  ifourier_ft(complex<double>* X, T* x, int N);

/****** Discrete sine transform ******/
template <typename T> void      sine_ft(double* X, T* x, int N);

/****** Thomas algorithm for solving linear systems ******/
void thomas_solver(double* x, double** M, double* b, int N);

#endif