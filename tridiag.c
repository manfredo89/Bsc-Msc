// parts of the function for backward Euler
void backward_euler(int xsteps, int tsteps, double delta_x, double alpha)
{
double *v, *r, a, b, c;
v = new double[xsteps+1]; // This is u
r = new double[xsteps+1]; // Right side of matrix equation Av=r
// Initialize vectors
for (int i = 0; i < xsteps; i++) {
r[i] = v[i] = func(delta_x*i);
}
r[xsteps] = v[xsteps] = 0;
// Matrix A, only constants
a = c = - alpha;
b = 1 + 2*alpha;
// Time iteration
for (int t = 1; t <= tsteps; t++) {
// here we solve the tridiagonal linear set of equations
tridag(a, b, c, r, v, x_steps+1);
// boundary conditions
v[0] = 0;
v[xsteps] = 0;
for (int i = 0; i <= x_steps; i++) {
r[i] = v[i];
}
}
...
}
// Function used to solve systems of equations for tridiagonal matrices
void tridag(double a, double b, double c, double *r, double *u, int n)
{
double bet, *gam;
gam = new double[n];
bet = b;
// forward substitution
u[0]=r[0]/bet;
for (int j=1;j<n;j++) {
gam[j] = c/bet;
bet = b - a*gam[j];
if (bet == 0.0) {cout << "Error 2 in tridag" << endl;}
u[j] = (r[j] - a*u[j-1])/bet;
}
// backward substitution
for (int j=n-2; j>=0; j--) {u[j] -= gam[j+1]*u[j+1];}
delete [] gam;
}
