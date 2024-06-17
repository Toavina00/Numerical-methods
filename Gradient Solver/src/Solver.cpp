#include "Solver.hpp"


template <typename T>
Solver<T>::Solver() {
    resolved = false;
}

template <typename T>
Solver<T>::Solver(int dim, T** M, T* b) {
    resolved = false;
    fill(dim, M, b);
}

template <typename T>
Solver<T>::Solver(string s) {
    resolved = false;
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
    delete[] r;

}

template <typename T>
void Solver<T>::fill(int dim, T** M, T* b) {

    this->dim = dim;
    this->M = M;
    this->b = b;

    x = new T[dim];
    y = new T[dim];
    r = new T[dim];

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
    r = new T[dim];

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
void Solver<T>::add(T* a, T* b, T* buff) {

    for (int i = 0; i < dim; i++) {
        buff[i] = a[i] + b[i];
    }

}

template <typename T>
void Solver<T>::sub(T* a, T* b, T* buff) {

    for (int i = 0; i < dim; i++) {
        buff[i] = a[i] - b[i];
    }

}

template <typename T>
void Solver<T>::scalar(T sc, T* a, T* buff) {

    for (int i = 0; i < dim; i++) {
        buff[i] = sc * a[i];
    }

}

template <typename T>
void Solver<T>::dot(T* a, T* buff) {
    
    for (int i = 0; i < dim; i++) {
        
        T tmp(0);
        
        for (int j = 0; j < dim; j++) {
            tmp += M[i][j] * a[j];
        }
    
        buff[i] = tmp;
    }

}


template <typename T>
T Solver<T>::product(T* a, T* b) {

    T res(0);

    for (int i = 0; i < dim; i++) {
        res += a[i]*b[i];
    }

    return res;

}

template <typename T>
T Solver<T>::norm(T* a) {
    
    return sqrt(product(a, a));

}

template <typename T> 
void Solver<T>::solve() {

    if (resolved) {
        return;
    }
    
    for (int i = 0; i < dim; i++) {
        x[i] = T(rand() % 100);
    }

    T eps(1/1000);

    T* tmp = new T[dim];
    
    dot(x, tmp);
    sub(b, tmp, r);

    while (norm(r) > eps) {
        

        dot(r, tmp);
        scalar((T) -1, tmp, y);

        T step = product(r, r) / product(y, r);

        scalar(step, r, tmp);
        sub(x, tmp, x);

        scalar(step, y, tmp);
        sub(r, tmp, r);

    }
    
    delete[] tmp;

    resolved = true;

}

template class Solver<float>;
template class Solver<double>;