#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <list>
#include <random>
#include <ctime>

using namespace std;

class CMK {

    /*
        ********************
        Cutthil-McKee solver
        ********************

        Let Mx = b a system of equation where
        M is symetric and sparse matrix, x and b are vectors.
        
        The algorithm solves the system as follow:
            (1) Find the Cuthill-McKee ordering that optimizes the profil
            (2) Permute the matrix and the vector b according to the ordering found
            (3) Store the profile into a vector AP
            (4) Compute the LDLt factorization according to the profil
            (5) Find Px the permuted solution by solving the equation LDLt(Px) = (Pb)
            (6) Retrieve the solution by re-ordering Px

    */

    private:
        vector<double> AP, LD, b, x;
        vector<vector<double>> M;
        vector<int> nDiag, order;
        int dim;


    public:
        /// Constructor and Destructor
        CMK()  {}
        ~CMK() {}

        /// fill the properties
        void fill (string path);

        /// Compute exentricity from a starting node
        int  eccentricity(int node);

        /// Build tree from a node
        list<list<int>> nodeTree(int node);

        /// Matrix setter and getter
        void set(int i, int j, double val);
        double get(int i, int j);

        /// Factorized Matrix setter and getter
        void dset(int i, int j, double val);
        double dget(int i, int j);

        /// Profile i_start
        int pstart(int i);

        /// Profile i_width
        int pwidth(int i);

        /// LDLt Factorisation
        void factorisation();

        /// Cuthill McKee ordering
        int cmkordering();

        /// Permutations
        void permuteMat(); // Matrix ordering
        void permuteVec(); // Vector ordering
        void permuteX();   // Solution re-ordering

        /// Store profil
        void storeProfil();

        /// Resolution
        void solve();

        /// Display
        void displayM();
        void displayA();
        void displayB();
        void displayL();
        void displayP();
        void solution();

};

template <typename T>
void display(vector<T> v) {
    int dim = v.size();
    for (int i(0); i < dim; i++) {
        cout << setw(8) << v[i];
    }
    cout << endl;
}

void CMK::displayM() {
    for (int i(0); i < dim; i++) {
        for (int j(0); j < dim; j++) {
            cout << setw(8) << M[i][j];
        }
        cout << endl;
    }
}

void CMK::displayA() {
    cout << setw(6) <<"AP:"; display(AP);  
    cout << setw(6) << "nDiag:"; display(nDiag);  
}

void CMK::displayL() {
    cout << setw(6) <<"Fact:"; display(LD);  
}

void CMK::displayP() {
    cout << setw(6) <<"Order:"; display(order);  
}

void CMK::displayB() {
    cout << setw(6) <<"b:"; display(b);  
}

void CMK::solution() {
    cout << setw(6) <<"x:"; display(x);  
}

int CMK::pwidth(int i) {
    if (i == 0) return 1;
    return nDiag[i] - nDiag[i-1];
}

int CMK::pstart(int i) {
    return i+1-pwidth(i);
}


void CMK::fill(string path) {

    fstream f(path, ios::in);

    if (!f.is_open()) {
        cerr << "Error: Cannot open the file" << endl;
        exit(1);
    }

    f >> dim;

    M = vector<vector<double>>(dim, vector<double>(dim));

    this->b = vector<double>(dim);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            f >> M[i][j];
        }
    }

    for (int i = 0; i < dim; i++) f >> b[i];
    
    f.close();

}

void CMK::set(int i, int j, double val) {
    if (i < j) return set(j, i, val);
    if (j < pstart(i)) return;
    AP[nDiag[i] - i + j] = val;
}

double CMK::get(int i, int j) {
    if (i < j) return get(j, i);
    if (j < pstart(i)) return 0;
    return AP[nDiag[i] - i + j];
}

void CMK::dset(int i, int j, double val) {
    if (i < j) return dset(j, i, val);
    if (j < pstart(i)) return;
    LD[nDiag[i] - i + j] = val;
}

double CMK::dget(int i, int j) {
    if (i < j) return dget(j, i);
    if (j < pstart(i)) return 0;
    return LD[nDiag[i] - i + j];
}

list<list<int>> CMK::nodeTree(int node) {

    vector<bool> visited(dim);
    vector<bool> peeked(dim);
    list<list<int>> N;
    list<int> temp;

    N.push_back(list<int>());
    N.back().push_back(node);

    list<int>::iterator it = N.back().begin();
    
    peeked[node] = true;
    visited[node] = true;

    int ecc = 0;

    while (it != N.back().end()) {
        
        int u = *it;

        for (int v(0); v < dim; v++) {
            if (M[u][v] && !visited[v] && !peeked[v]) {
                peeked[v] = true;
                temp.push_back(v);
            }
        }

        visited[u] = true;
        it++;

        if (it == N.back().end()) {
            if (temp.empty()) break;
            else {
                N.push_back(list<int>());
                while (!temp.empty()) {
                    N.back().push_back(temp.back());
                    temp.pop_back();
                }
                it = N.back().begin();
            }
        }
    }

    return N;
}

int CMK::eccentricity(int node) {
    return nodeTree(node).size() - 1;
}

void CMK::factorisation() {

    LD = vector<double>(AP.size());

    for (int i(0); i < dim; i++) {
        
        double sum(0.0); 
        for (int k = pstart(i); k < i; k++) 
            sum += dget(k, k) * dget(i, k) * dget(i, k);
        dset(i, i, get(i, i) - sum); 
        
        for (int j = pstart(i); j < i; j++) {
            double s(0.0);
            for (int k = pstart(j); k < j; k++)
                s += dget(i, k) * dget(k, k) * dget(j, k);
            dset(i, j, (get(i, j)-s)/dget(j, j));
        }
    }
    
}

int CMK::cmkordering() {

    srand(time(NULL));
    int n = rand() % dim;

    vector<bool> visited(dim);
    list<list<int>> nTree;
    bool end = false;
    
    while (!end) {
        
        visited[n] = true;
        nTree = nodeTree(n);
        
        end = true;

        for (auto s: nTree.back()) {
            if (!visited[s]) {
                list<list<int>> sTree = nodeTree(s);
                if (nTree.size() < sTree.size()) {
                    n = s;
                    end = false;
                    break;    
                }
            }
        }

    }

    order.clear();

    for (auto& level: nTree) {
        //level.sort([this] (int a, int b) {
        //    int sa(0), sb(0);
        //    for (int i(0); i < this->dim; i++) {
        //        if (this->M[a][i]) sa++;
        //        if (this->M[b][i]) sb++;
        //    }
        //    return sa < sb;
        //});
        order.insert(order.end(), level.begin(), level.end());
    }

    for (int i(0); i < dim / 2 ; i++) swap(order[i], order[dim - i - 1]);

    return n;
}

void CMK::storeProfil() {
    
    this->AP = vector<double>();
    this->nDiag = vector<int>(dim);

    for (int i(0); i < dim; i++) {
        bool skip = true;
        for (int j(0); j <= i; j++) {
            if (M[i][j] != 0) skip = false;
            if (skip) continue; 
            AP.push_back(M[i][j]);
        }
        nDiag[i] = AP.size() - 1;
    }

}

void CMK::permuteMat() {
    vector<vector<double>> Tmp(dim);
    for (int i(0); i < dim; i++) {
        Tmp[i] = M[order[i]];
    }
    for (int i(0); i < dim; i++) {
        for (int k(0); k < dim; k++) {
            M[i][k] = Tmp[i][order[k]];
        }
    }
}

void CMK::permuteVec() {
    vector<double> tmp(dim);
    for (int i(0); i < dim; i++) tmp[i] = b[order[i]];
    b = tmp;
}   

void CMK::permuteX() {
    vector<double> tmp(dim);
    vector<int> invOrder(dim);
    for (int i(0); i < dim; i++) invOrder[order[i]] = i;
    for (int i(0); i < dim; i++) tmp[i] = x[invOrder[i]];
    x = tmp;
}   

void CMK::solve() {

    x = vector<double>(dim, 0.0);

    for (int k = 0; k < dim; k++) {
        x[k] = b[k];
        for (int i = pstart(k); i < k; i++) {
            x[k] -= dget(k, i) * x[i];
        }
    }
    
    for (int k = 0; k < dim; k++) x[k] /= dget(k, k);

    for (int k = dim - 1; k >= 0; k--) {
        for (int i = k+1; i < dim; i++) {
            x[k] -= dget(i, k) * x[i];
        }
    }
}

int main() {

/// Initialisation
    CMK mat;
    string path("syst.txt");

    cout << "----------------------" << endl;
    cout << " Cuthill-Mckee solver " << endl;
    cout << "----------------------" << endl;

    cout << endl;


/// Solve

    // Input the system into the program
    mat.fill(path);
    
    fixed(cout);
    cout << setprecision(2);

    cout << "----Display matrix----" << endl;
    mat.displayM();
    cout << "----------------------" << endl;

    cout << "----Display vector----" << endl;
    mat.displayB();
    cout << "----------------------" << endl;

    cout << endl;


    // Find the Cuthill-McKee ordering
    int n = mat.cmkordering();
    list<list<int>> levels = mat.nodeTree(n);

    // Permute the matrix and the vector b
    mat.permuteMat();
    mat.permuteVec();

    // Store the profile
    mat.storeProfil();

    // Compute the LDLt factorization
    mat.factorisation();

    // Find the permuted solution
    mat.solve();

    // Re-order (Px) to get the solution
    mat.permuteX();


/// Displays

    cout << "---------Start--------" << endl;
    cout << "CMK Start: " << n << endl;
    cout << "Eccentricity: " << levels.size() - 1 << endl;
    cout << "----------------------" << endl;
    cout << "---------Levels--------" << endl;
    
    list<list<int>>::iterator it0;
    list<int>::iterator it1;

    for (it0 = levels.begin(); it0 != levels.end(); it0++) {
        cout << "[ ";
        for (it1 = it0->begin(); it1 != it0->end(); it1++)  
            cout << *it1 << " ";
        cout << "]" << endl;
    }

    cout << "----------------------" << endl; 

    cout << "-----Display order----" << endl;
    mat.displayP();
    cout << "----------------------" << endl;

    cout << "----Display profil----" << endl;
    mat.displayA();
    cout << "----------------------" << endl;

    cout << "--Display factorized--" << endl;
    mat.displayL();
    cout << "----------------------" << endl;

    cout << "-------Solution-------" << endl;
    mat.solution();
    cout << "----------------------" << endl;

    return 0; 
}