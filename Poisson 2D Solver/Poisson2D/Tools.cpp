#include "Tools.hpp"

void display(double* x, int N) {
    fixed(cout);
    cout << setprecision(4);
    for (int i(0); i < N; i++) {
        cout << setw(10) << x[i];
    }
    cout << endl;
}

void display(complex<double>* x, int N) {
    fixed(cout);
    cout << setprecision(4);
    for (int i(0); i < N; i++) {
        cout << setw(10) << x[i].real() << " + j" << x[i].imag();
    }
    cout << endl;
}

void display(double** x, int N, int M) {
    fixed(cout);
    cout << setprecision(4);
    for (int i(0); i < N; i++) {
        for (int j(0); j < M; j++) {
            cout << setw(10) << x[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

template <typename T> void computation(complex<double> *y, T *x, int N, int s, bool forward) {
    if (N == 1) {
    /// Return the unique term if trying to compute the fourier transformation of an array of size one
        y[0] = x[0];   
    } else {
    
    /// Split into even index and odd index sub arrays
        computation(y, x, N/2, 2*s, forward);
        computation(&y[N/2], &x[s], N/2, 2*s, forward);
    
    /// Compute fourier transformation of sub arrays
        for (int i(0); i < N / 2; i++) {
            complex<double> w = (forward)? 
                complex<double> (cos(-2*M_PI*i/N), sin(-2*M_PI*i/N)):
                complex<double> (cos( 2*M_PI*i/N), sin( 2*M_PI*i/N));
            complex<double> p = y[i];
            complex<double> q = y[i + N/2] * w;
            y[i] = p + q;
            y[i + N/2] = p - q;
        }
    }
}

template <typename T> void factor(T* X, double x, int N) {
    for (int i(0); i < N; i++) {
        X[i] *= x;
    }
}

template <typename T> void fourier_ft(complex<double>* X, T* x, int N) {
    computation(X, x, N, 1, true);
}

template <typename T> void ifourier_ft(complex<double>* X, T* x, int N) {
    computation(X, x, N, 1, false);
    factor(X, 1.0/N, N);
}

template <typename T> void sine_ft(double* X, T* x, int N) {
/// Initialisation    
    T*  tmp = new T [N];
    complex<double>* Y = new complex<double>[N];
    double p = sqrt((double) N/2); 

/// Compute the auxiliary array
    tmp[0] = 0;
    for (int i(1); i < N; i++) {
        tmp[i] = (x[i]-x[N-i])*0.5 + (x[i]+x[N-i])*sin(i*M_PI/N);
    }

/// Compute the fast fourier transformation inverse of the auxiliary array
    ifourier_ft(Y, tmp, N);

/// Get the discrete sinus transformation
    X[1] = p * Y[0].real();
    for (int i(1); i < N/2; i++) {
        X[2*i] = 2 * p * Y[i].imag();
        X[2*i+1] = 2 * p * Y[i].real() + X[2*i-1];
    }

/// Free memory
    delete [] tmp;
    delete [] Y;
}

void thomas_solver(double *x, double** M, double* b, int dim) {
/// Initialisation 
    double *c, *d;
    c = new double[dim];    
    d = new double[dim];    

/// Upper Triangularisation 
    for (int i(0); i < dim; i++) {
        if (i == 0) {
            c[i] = M[i][i+1] / M[i][i];
            d[i] = b[i] / M[i][i];
        } else if (i == (dim - 1)) {
            d[i] = (b[i] - M[i][i-1] * d[i-1]) / (M[i][i] - M[i][i-1]*c[i-1]);
        } else {
            c[i] = M[i][i+1] / (M[i][i] - M[i][i-1]*c[i-1]);
            d[i] = (b[i] - M[i][i-1] * d[i-1]) / (M[i][i] - M[i][i-1]*c[i-1]);
        }
    }


/// Resolution of the system
    for (int k = dim - 1; k >= 0; k--) {
        if (k == dim - 1) {
            x[k] = d[k];
        } else {
            x[k] = d[k] - c[k] * x[k + 1];
        }
    }
    
/// Free allocated memory
    delete[] c;
    delete[] d;
}


template void       factor(double* X, double x, int N);
template void      sine_ft(double* X, double* x, int N);
template void   fourier_ft(complex<double>* X, double* x, int N);
template void  ifourier_ft(complex<double>* X, double* x, int N);

template void       factor(complex<double>* X, double x, int N);
template void      sine_ft(double* X, complex<double>* x, int N);
template void   fourier_ft(complex<double>* X, complex<double>* x, int N);
template void  ifourier_ft(complex<double>* X, complex<double>* x, int N);