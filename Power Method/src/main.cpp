#include "PowerMethod.hpp"


int main() {

    cout << "System of equation solver" << endl 
        << "---------------------------" << endl;

/// Variables' declarations
    string path;
    PowerMethod P;

/// Initialisation
    cout << "Give the path to the file: ";
    getline(cin, path);
    P.fill(path);

/// Treatement and display
    P.display();
    P.solution();

    float* a = new float[2];
    float* b = new float[2];

    a[0] = 10;
    b[0] = 20;

    P.swap(&a, &b);
    cout << a[0] << " " << b[0] << endl;
}