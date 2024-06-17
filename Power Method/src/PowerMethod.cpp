#include "PowerMethod.hpp"


PowerMethod::PowerMethod() {
}

PowerMethod::PowerMethod(int dim, float** M) {
    fill(dim, M);
}

PowerMethod::PowerMethod(string s) {
    fill(s);
}

PowerMethod::~PowerMethod() {

    for(int i = 0; i < dim; i++) {
        delete[] M[i];
        delete[] U[i];
    }
    delete[] M;
    delete[] x;

}

void PowerMethod::fill(int dim, float** M) {

    this->dim = dim;
    this->M = M;

    x = new float[dim];

    U = new float*[dim];
    for (int i = 0; i < dim; i++) {
        U[i] = new float[dim];
    }
}

void PowerMethod::fill(string path) {

    fstream f(path, ios::in);

    if (!f.is_open()) {
        cerr << "Error: Cannot open the file" << endl;
        exit(1);
    }

    f >> dim;

    M = new float*[dim];
    x = new float[dim];

    for (int i = 0; i < dim; i++) {
        M[i] = new float[dim];
        for (int j = 0; j < dim; j++) {
            f >> M[i][j];
        }
    }

    U = new float*[dim];
    for (int i = 0; i < dim; i++) {
        U[i] = new float[dim];
    }
    
    f.close();

}

void PowerMethod::display(int prec) {

    cout << "The Matrix:" << endl;
    
    cout << fixed << setprecision(prec);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            cout << setw(8+prec) << M[i][j];
        }
        cout << endl;
    }

}

void PowerMethod::solution(int prec) {

    cout << "The eigenvalues:" << endl;

    cout << fixed << setprecision(prec);

    compute(0);

    for (int i = 0; i < dim; i++) {
        cout << setw(8+prec) << x[i];
        break;
    }

    cout << endl;

    cout << "The eigenvectors:" << endl;
    
    for (int i = 0; i < dim; i++) {
        cout << setw(8+prec) << U[0][i];
    }
    cout << endl;

}

void PowerMethod::dot(float* v, float* b) {
    for (int i = 0; i < dim; i++) {
        b[i] = 0;
        for (int j = 0; j < dim; j++) {
            b[i] += M[i][j] * v[j];
        }
    }
}


void PowerMethod::normalize(float* v) {
    float nrm = norm(v);
    for (int i = 0; i < dim; i++) {
        v[i] /= nrm;
    }
}


void PowerMethod::swap(float** v1, float** v2) {
    float *tmp;
    tmp = *v1;
    *v1 = *v2;
    *v2 = tmp;
}


float PowerMethod::product(float* v1, float* v2) {
    float s(0);
    for (int i = 0; i < dim; i++) {
        s += v1[i]*v2[i];
    }
    return s;
}


float PowerMethod::norm(float* v) {
    return sqrt(product(v, v));
}


void PowerMethod::compute(int k) {
    float a, b;
    float *u, *v;
    
    u = new float[dim];
    v = new float[dim];

    for (int i = 0; i < dim; i++) {
        u[i] = (float) (rand() % 10);
    }

    normalize(u);
    dot(u, v);

    a = product(u, v);
    
    swap(&u, &v);
    normalize(u);
    dot(u, v);

    b = product(u, v);

    float eps(1/100000);

    while (abs(b - a) > eps) {
        
        a = b;
        swap(&u, &v);
        
        normalize(u);
        dot(u, v);

        b = product(u, v);
    
    }

    x[k] = a;
    U[k] = u;

    delete[] v;
    
}