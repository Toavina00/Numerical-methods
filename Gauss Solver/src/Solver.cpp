#include "Solver.hpp"


template <typename T>
Solver<T>::Solver() {
    resolved = false;
    triangular = false;
}

template <typename T>
Solver<T>::Solver(int dim, T** M, T* b) {
    resolved = false;
    triangular = false;
    fill(dim, M, b);
}

template <typename T>
Solver<T>::Solver(string s) {
    resolved = false;
    triangular = false;
    fill(s);
}

template <typename T>
Solver<T>::~Solver() {

    for(int i = 0; i < dim; i++) {
        delete[] M[i];
    }
    delete[] M;
    delete[] b;
    delete[] x;

}

template <typename T>
void Solver<T>::fill(int dim, T** M, T* b) {

    this->dim = dim;
    this->M = M;
    this->b = b;

    x = new T[dim];

}

template <typename T>
void Solver<T>::fill(string path) {

    fstream f(path, ios::in);

    if (!f.is_open()) {
        cerr << "Error: Cannot open the file" << endl;
        exit(1);
    }

    f >> dim;

    M = new T*[dim];
    b = new T[dim];
    x = new T[dim];

    for (int i = 0; i < dim; i++) {
        M[i] = new T[dim];
        for (int j = 0; j < dim; j++) {
            f >> M[i][j];
        }
        char del; f >> del;
        f >> b[i];
    }
    
    f.close();

}

template <typename T>
void Solver<T>::display(int prec) {

    cout << "The system of equation:" << endl;
    
    cout << fixed << setprecision(prec);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            cout << setw(8+prec) << M[i][j];
        }
        cout << setw(8+prec) << "|";
        cout << setw(8+prec) << b[i] << endl;
    }

}

template <typename T>
void Solver<T>::solution(int prec) {

    cout << "The solution to ths system of equation:" << endl;

    if (!resolved) {
        solve();
    }

    cout << fixed << setprecision(prec);

    for (int i = 0; i < dim; i++) {
        cout << setw(8+prec) << "x" << i << setw(8+prec) << "=" << setw(8+prec) << x[i] << endl;
    }

}

template <typename T>
void Solver<T>::swap(int r0, int r1) {
    
    T* tmp = M[r0];
    M[r0] = M[r1];
    M[r1] = tmp;

    T t = b[r0];
    b[r0] = b[r1];
    b[r1] = t;

}

template <typename T>
void Solver<T>::upperTri() {

    if (triangular) {
        return;
    }

    for (int k = 0; k < dim - 1; k++) {

        int pivot = k;
        for (int i = k; i < dim; i++) {
            if (abs(M[pivot][k]) < abs(M[i][k])) {
                pivot = i;
            }
        }
        if ( pivot != k) {
            swap(k, pivot);
        }

        for (int i = k+1; i < dim; i++) {
            for (int j = k+1; j < dim; j++) {
                M[i][j] = M[i][j] - (M[i][k]/M[k][k])*M[k][j];
            }
            b[i] = b[i] - (M[i][k]/M[k][k])*b[k];
            M[i][k] = 0;
        }

    }

    triangular = true;

}

template <typename T> 
void Solver<T>::solve() {

    if (resolved) {
        return;
    }
    
    if (!triangular) {
        upperTri();
    }

    for (int k = dim - 1; k >= 0; k--) {
        T s(0);
        for (int i = k+1; i < dim; i++) {
            s += M[k][i] * x[i];
        }
        x[k] = b[k] - s;
        x[k] /= M[k][k];
    }

    resolved = true;

}

template class Solver<float>;
template class Solver<double>;