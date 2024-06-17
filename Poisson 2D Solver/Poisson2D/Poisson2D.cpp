#include "Poisson2D.hpp"

Poisson2D::Poisson2D() {
}

Poisson2D::~Poisson2D() {
}

void Poisson2D::setDomain(double x, double x_len, double y, double y_len, int M, int N) {
    this->x_len = x_len;
    this->y_len = y_len;
    this->x = x;
    this->y = y;
    this->M = M;
    this->N = N;
    this->dy = y_len/(N+1);
    this->dx = x_len/(M+1);
}


function<double(int, int)> Poisson2D::discretize(function<double(double, double)> func) {
    return [this, func] (int i, int j) {
        return func(this->x + i*this->dx, this->y + j*this->dy);
    };
}


double** Poisson2D::solve(function<double(double, double)> f_func, function<double(double, double)> g_func) {
/// Initialisation
    double *d;
    double a, b, c;
    double **n, **u, **s, **h, **v, **T, **y;

/// Function discretization
    function<double(int, int)> f = discretize(f_func);
    function<double(int, int)> g = [g_func, this] (int i, int j) {
        if (i == 0 || i == (this->M+1) || j == 0 || j == (this->N+1)) return discretize(g_func)(i,j);
        else return 0.0;
    };

/// Computing constants
    a = -(dy*dy)/(dx*dx);
    b = (dy*dy)*(2/(dx*dx)+2/(dy*dy));
    c = -1;

/// Computing the eigenvalues
    d = new double[M];
    for (int i(0); i < M; i++) {
        d[i] = 2 + 2 * (dy/dx) * (dy/dx) * (1 - cos((i+1)*M_PI/(M+1)));
    }


/// Computing S_ij
    s = new double*[N];
    s--;
    for (int i(1); i <= N; i++) {
        s[i] = new double[M];
        s[i]--;
        for (int j(1); j <= M; j++) {
            s[i][j] = (dy*dy)*f(j, i) - c*(g(j-1, i) + g(j+1, i)) - a*(g(j, i-1) + g(j, i+1));
        }
        s[i]++;
    }
    s++;


/// Computing h_j
    n = new double*[N];
    for (int i(0); i < N; i++) {
        n[i] = new double[M];
        sine_ft(n[i] - 1, s[i] - 1, M + 1);
    }

    h = new double*[M];
    for (int j(0); j < M; j++) {
        h[j] = new double[N];
        for (int i(0); i < N; i++) {
            h[j][i] = n[i][j];
        }
    }


/// Solving the equation T_j @ u_j = h_j
    T = new double*[N];
    for (int i(0); i < N; i++) {
        T[i] = new double[N];
        for (int j(0); j < N; j++) {
            T[i][j] = 0;
        }
    }
    
    u = new double*[M];
    for (int j(0); j < M; j++) {
        u[j] = new double[N];
        for (int i(0); i < N; i++) {
            T[i][i] = d[j];
            if (i-1 >= 0) {
                T[i][i-1] = -1;
            }
            if (i+1 < N) {
                T[i][i+1] = -1;
            }
        }
        thomas_solver(u[j], T, h[j], N);
    }


/// Computing v_i
    v = new double*[N];
    for (int i(0); i < N; i++) {
        v[i] = new double[M];
        for (int j(0); j < M; j++) {
            v[i][j] = u[j][i];
        }
    }


/// Computing y_i
    y = new double*[N];
    for (int i(0); i < N; i++) {
        y[i] = new double[M];
        sine_ft(y[i] - 1, v[i] - 1, M + 1);
    }
    

/// Freeing memory
    for (int i(0); i < N; i++) {
        delete[] n[i];
        delete[] v[i];
        delete[] s[i];
        delete[] T[i];
    }

    for (int j(0); j < M; j++) {
        delete[] h[j];
        delete[] u[j];
    }

    delete[] n;
    delete[] v;
    delete[] s;
    delete[] T;
    delete[] h;
    delete[] u;
    delete[] d;


    return y;
}