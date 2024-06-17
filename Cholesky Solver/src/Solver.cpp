#include "Solver.hpp"


template <typename T>
Solver<T>::Solver() {
    resolved = false;
    factored = false;
}

template <typename T>
Solver<T>::Solver(int dim, T** M, T* b) {
    resolved = false;
    factored = false;
    fill(dim, M, b);
}

template <typename T>
Solver<T>::Solver(string s) {
    resolved = false;
    factored = false;
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
    delete[] y;

}

template <typename T>
void Solver<T>::fill(int dim, T** M, T* b) {

    this->dim = dim;
    this->M = M;
    this->b = b;

    x = new T[dim];
    y = new T[dim];

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
    y = new T[dim];

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

    cout << "The matrix:" << endl;
    
    cout << fixed << setprecision(prec);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            cout << setw(8+prec) << M[i][j];
        }
        cout << setw(8+prec) << "|";
        cout << setw(8+prec) << b[i];
        cout << endl;
    }

}

template <typename T>
void Solver<T>::solution(int prec) {

    cout << "The solution to the system of equation:" << endl;

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
void Solver<T>::factorization() {

    if (factored) {
        return;
    }

    for (int i = 0; i < dim; i++) {

        T s(0);
        
        for (int j = 0; j < i; j++) {

            s = 0;
            for (int k = 0; k < j; k++) {
                s += M[i][k] * M[j][k];
            }

            M[i][j] -= s;
            M[i][j] /= M[j][j];

        }

        s = 0;
        for (int j = 0; j < i; j++) {
            s += pow(M[i][j], 2);          
        }
        
        M[i][i] -= s;
        M[i][i] = pow(M[i][i], 0.5);

        //for (int j = i + 1; j < dim; j++) {
        //    M[i][j] = 0;
        //}

    }

    factored = true;

}

template <typename T> 
void Solver<T>::solve() {

    if (resolved) {
        return;
    }
    
    if (!factored) {
        factorization();
    }

    for (int k = 0; k < dim; k++) {
        T s(0);
        for (int i = 0; i < k; i++) {
            s += M[k][i] * x[i];
        }
        x[k] = b[k] - s;
        x[k] /= M[k][k];
    }

    for (int k = dim - 1; k >= 0; k--) {
        T s(0);
        for (int i = k+1; i < dim; i++) {
            s += M[i][k] * x[i];
        }
        x[k] = x[k] - s;
        x[k] /= M[k][k];
    }

    resolved = true;

}

template class Solver<float>;
template class Solver<double>;