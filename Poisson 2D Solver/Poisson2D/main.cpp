#include <streambuf>

#include "Tools.hpp"
#include "Poisson2D.hpp"

double f_func(double x, double y);
double g_func(double x, double y);
double u_func(double x, double y);

int main() {
    // Initialisation
    Poisson2D solver;
    double **y, **ty, **err, error;


    // Set up the grid and boundary conditions
    int N;                    // Discretisation grid along the x-axis
    int M;                    // Discretisation grid along the y-axis
    double a, b, c, d;        // Discretisation domain

    // Initialize the solver
    cout << "---------------------------" << endl;
    cout << "---- Poisson 2d solver ----" << endl;
    cout << "---------------------------" << endl;
    cout << endl;
    cout << "Let the following system:" << endl;
    cout << "\t -u\"(x,y) = -5.0 * (x*x + y*y - 2) * exp(-(0.5)*(x*x+y*y))" << endl;
    cout << endl;
    cout << "We will solve the system in D = [a, b] x [c, d] on M x N points" << endl;
    cout << endl;
    cout << "Enter a: "; cin >> a;
    cout << "Enter b: "; cin >> b;
    cout << "Enter c: "; cin >> c;
    cout << "Enter d: "; cin >> d;
    cout << "Enter the number of points: " << endl;
    cout << "Enter M: "; cin >> M;
    cout << "Enter N: "; cin >> N;
    
    double sx(a), sy(c), Lx(b-a), Ly(d-c);
    solver.setDomain(sx, Lx, sy, Ly, M, N);

    // Discretize function
    function<double(int, int)> u = solver.discretize(u_func);


    // Solve the problem
    y = solver.solve(f_func, g_func);


    // Printing results
    fixed(cout);
    cout << setprecision(2);

    cout << "---------------" << endl;
    cout << "Output" << endl;
    cout << "---------------" << endl;
    display(y, N, M);

    cout << "---------------" << endl;
    cout << "Expected Output" << endl;
    cout << "---------------" << endl;
    ty  = new double* [N];
    err = new double* [N];
    for (int i(0); i < N; i++) {
        ty[i]  = new double[M];
        err[i] = new double[M];
        for (int j(0); j < M; j++) {
            ty[i][j]  = u(j+1, i+1);
            err[i][j] = abs(ty[i][j]-y[i][j]);
            error += (ty[i][j]-y[i][j])*(ty[i][j]-y[i][j]);
        }
    }
    display(ty, N, M);

    cout << "---------------" << endl;
    cout << "Absolute error" << endl;
    cout << "---------------" << endl;
    display(err, N, M);


    cout << "---------------" << endl;
    cout << "RMS Error: " << sqrt(error) << endl;
    cout << "---------------" << endl;

    /*
        The following code is used to plot the results
        Require Gnuplot
    

    // Plot the results
    FILE* gnuplot = popen("gnuplot -persist", "w");
    
    if  (!gnuplot) {
        cerr << "Error opening gnuplot!" << endl;
        exit(EXIT_FAILURE);
    }

    stringstream gp;

    gp << "set term wxt size 640,480\n";
    gp << "set term wxt title 'Solution'\n";
    gp << "set nokey\n";
    gp << "set xzeroaxis\n";
    gp << "set yzeroaxis\n";
    gp << "set autoscale z\n";
    gp << "set yrange ["<<c<<":"<<d<<"]\n";
    gp << "set xrange ["<<a<<":"<<b<<"]\n";

    gp << "$points << EOD\n";
    double dx(Lx/(M+1)), dy(Ly/(N+1));
    for (int i(0); i < N; ++i) {
        for (int j(0); j < M; ++j) {
            gp  << sx + (i+1)*dx << " " 
                << sy + (j+1)*dy << " "
                << y[i][j]   << " "
                << ty[i][j]  << " "
                << err[i][j] << "\n";
        }
    }
    gp << "EOD\n";

    gp << "set xlabel 'X'\n"
       << "set ylabel 'Y'\n"
       << "set zlabel 'Z(x,y)'\n";

    gp << "splot" << " ";
    gp << "$points u 1:2:3 with pm3d title 'approximated u(x, y)'" << " , ";
    gp << "$points u 1:2:4 with pm3d title 'true u(x, y)'" << " , ";
    gp << "$points u 1:2:5 with pm3d title 'error approximation'" << "\n";

    fprintf(gnuplot, "%s\n", gp.str().c_str());

    fflush(gnuplot);
    pclose(gnuplot);

    // Free memory
    for (int i(0); i < N; i++) {
        delete []  ty[i];
        delete []   y[i];
        delete [] err[i];
    }
    delete[]  ty;
    delete[]   y;
    delete[] err;

    */

    return 0;
}

double f_func(double x, double y){
    // Here is defined the f function
    return -5.0 * (x*x + y*y - 2) * exp(-(0.5)*(x*x+y*y));
}

double u_func(double x, double y){
    // Here is defined the expected u function
    return 5.0 * exp(-(0.5)*(x*x+y*y));
}

double g_func(double x, double y) {
    // Here is defined the g function
    return 5.0 * exp(-(0.5)*(x*x+y*y));
}