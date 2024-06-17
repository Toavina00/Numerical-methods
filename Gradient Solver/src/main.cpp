#include "Solver.hpp"


int main() {

    cout << "Gradient Descend" << endl 
        << "---------------------------" << endl;

    srand(time(nullptr));

/// Variables' declarations
    string path;
    Solver<double> S;

/// Initialisation
    cout << "Give the path to the file: ";
    getline(cin, path);
    S.fill(path);

/// Treatement and display
    S.display();
    S.solution();
}