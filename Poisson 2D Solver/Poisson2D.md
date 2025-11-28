# 2D Poisson Problem

## Presentation of the 2D Poisson Problem
The 2D Poisson problem is an elliptic partial differential equation that appears in many fields of physics and engineering. It involves finding a function $u$ defined over a domain $\Omega \in \mathbb{R}^2$ that satisfies the following equation:

$$ 
\begin{cases}
-\Delta u(x, y) = f(x, y) \\
u(x, y) = g(x, y) \quad \text{on} \quad \partial \Omega 
\end{cases}
$$

## Solving the 2D Poisson Problem
To solve the problem, we need to discretize the study domain $\Omega$ into a grid of $(M+2) \times (N+2)$ regularly spaced points. Let $dx$ and $dy$ be the discretization steps in $x$ and $y$, respectively. Knowing the values of $u(x, y)$ on $\partial \Omega$, we propose to find approximations of $u(x, y)$ on the interior points of the domain $(x_i, y_j)_{1 \leq i \leq M, 1 \leq j \leq N}$. By developing our system, we obtain:

$$
\begin{cases}
-\frac{1}{dx^2}(u_{i,j-1}-2u_{i,j}+u_{i,j+1}) - \frac{1}{dy^2}(u_{i-1,j}-2u_{i,j}+u_{i+1,j}) = f_{i,j} \\
u_{i,j} = g_{i,j} \quad \text{if} \quad i = 0 \ \text{or} \ i = N+1 \ \text{or} \ j = 0 \ \text{or} \ j = M+1 \\
0 \leq i \leq N+1 \quad \text{and} \quad 0 \leq j \leq M+1 \\
\end{cases}
$$

**Notation**: For any function $f$, we have $f_{i,j} = f(x_j, y_i)$. Defining $\alpha = -\frac{1}{dx^2}$, $\beta = \frac{2}{dx^2} + \frac{2}{dy^2}$, $\gamma = -\frac{1}{dy^2}$, the equation can be generally written as:

$$
\gamma u_{i-1,j} + \alpha u_{i,j-1} + \beta u_{i,j} + \alpha u_{i,j+1} + \gamma u_{i+1,j} = f_{i,j}
$$

Given that $u_{i,j} = g_{i,j}$ if $i = 0$ or $i = N+1$ or $j = 0$ or $j = M+1$, we obtain the following system:

$$
\text{For} \ i=1 \begin{cases}
\beta u_{1,1} + \alpha u_{1,2} + \gamma u_{2,1} &= f_{1,1} - \gamma g_{0,1} - \alpha g_{1,0} \quad \text{for} \ j=1 \\
\alpha u_{1,j-1} + \beta u_{1,j} + \alpha u_{1,j+1} + \gamma u_{2,j} &= f_{1,j} - \gamma g_{0,j} \quad \text{for} \ j=2,...,M-1\\
\alpha u_{1,M-1} + \beta u_{1,M} + \gamma u_{2,M} &= f_{1,M} - \gamma g_{0,M} - \alpha g_{1,M+1} \quad \text{for} \ j=M\\
\end{cases} 
$$
$$
\text{For} \ i=2,...,N-1 \begin{cases}
\gamma u_{i-1,1} + \beta u_{i,1} + \alpha u_{i,2} + \gamma u_{i+1,1} &= f_{i,1} - \alpha g_{i,0} \quad \text{for} \ j=1\\
\gamma u_{i-1,j} + \alpha u_{i,j-1} + \beta u_{i,j} + \alpha u_{i,j+1} + \gamma u_{i+1,j} &= f_{i,j} \quad \text{for} \ j=2,...,M-1\\
\gamma u_{i-1,M} + \alpha u_{i,M-1} + \beta u_{i,M} + \gamma u_{i+1,M} &= f_{i,M} - \alpha g_{i,M+1} \quad \text{for} \ j=M\\
\end{cases} 
$$

$$
\text{For} \ i=N \begin{cases}
\gamma u_{N-1,1} + \beta u_{N,1} + \alpha u_{N,2} &= f_{N,1} - \alpha g_{N,0} - \gamma g_{N+1,1} \quad \text{for} \ j=1\\
\gamma u_{N-1,j} + \alpha u_{N,j-1} + \beta u_{N,j} + \alpha u_{N,j+1} &= f_{N,j} - \gamma g_{N+1,j} \quad \text{for} \ j=2,...,M-1\\
\gamma u_{N-1,M} + \alpha u_{N,M-1} + \beta u_{N,M} &= f_{N,M} - \alpha g_{N,M+1} - \gamma g_{N+1,M} \quad \text{for} \ j=M\\
\end{cases} 
$$

This can be written as:

$$ 
\begin{pmatrix}
B & C &   &   &   \\
C & .. &   &   &   \\
  &   & .. &   &   \\
  &   &   & .. & C \\
  &   &   & C & B \\
\end{pmatrix}
\begin{pmatrix}
Y_1 \\ Y_2 \\ .. \\ .. \\ Y_{N}
\end{pmatrix} =  
\begin{pmatrix}
S_1 \\ S_2 \\ .. \\ .. \\ S_{N}
\end{pmatrix}
$$

$$
\text{with} \quad B =
\begin{pmatrix}
\beta & \alpha &   &   \\
\alpha &  .. &   &   \\
  &   & .. &   &   \\
  &   & & ..  & \alpha \\
  &   & & \alpha & \beta \\
\end{pmatrix} \quad  C = \gamma I_m \quad \text{and} \quad S_i = \begin{pmatrix} S_{i,1} \\ S_{i,2} \\ .. \\ .. \\ S_{i,M} \end{pmatrix} 
$$

To solve the system, we propose using the **Fourier Tridiagonal Method**.

## Fourier Tridiagonal Method
Our system can be written as:

$$
\begin{cases}
            BY_1 + CY_2     &= S_1 \\
   CY_{i-1} + BY_i + CY_{i+1} &= S_i \quad 2 \leq i \leq N-1\\
   CY_{N-1} + BY_N           &= S_N \\
\end{cases}
$$

We know that $B$ is symmetric, so it is diagonalizable, and if we denote $Q$ as the change of basis matrix, we get:

$$
\begin{cases}
    D = Q^T B Q \\
    \Delta = Q^T C Q
\end{cases} \Rightarrow
\begin{cases}
    D Q^T = Q^T B \\
    \Delta Q^T = Q^T C
\end{cases}
$$

Let $\{d_j\}_{1 \leq j \leq M}$ and $\{\delta_j\}_{1 \leq j \leq M}$ be the eigenvalues of $D$ and $\Delta$ respectively, and by noting $Q^T S_i = \eta_i$ and $Q^T Y_i = u_i$, we get:

$$
\begin{cases}
            Du_1 + \Delta u_2     &= \eta_1 \\
   \Delta u_{i-1} + Du_i + \Delta u_{i+1} &= \eta_i \quad 2 \leq i \leq N-1\\
   \Delta u_{N-1} + Du_N           &= \eta_N \\
\end{cases}
$$

Then, by grouping our system by rows of the same index, we can write:

$$
\begin{cases}
    d_j u_{1,j} + \delta_j u_{2,j} &= \eta_{1,j} \\
    \delta_j u_{i-1,j} + d_j u_{i,j} + \delta_j u_{i+1,j} &= \eta_{i,j} \quad 2 \leq i \leq N-1 \\
    \delta_j u_{N-1,j} + d_j u_{N,j} &= \eta_{N,j} \\
\end{cases}_{1 \leq j \leq M}
$$

This corresponds to solving $M$ tridiagonal linear systems of the following form:

$$
T_j \nu_j = h_j \quad \quad j=1,...,M
$$

$$
\text{where} \quad T_j = 
\begin{pmatrix}
d_j & \delta_j & & & \\
\delta_j & & & & \\
& & & & \\
& & & & \delta_j \\
& & & \delta_j & d_j \\
\end{pmatrix} \quad \text{and} \quad \nu = u^t, \quad \eta = h^t
$$

Once the $M$ systems are solved, we obtain the solution to the problem as follows:

$$
\begin{cases}
Y_i = Qu_i \\
i = 1,...,N
\end{cases}
$$

Those systems are tridiagonal, which makes solving them relatively simple using the Thomas algorithm.

## Algorithm

Before solving the system, we first multiply it by $dy^2$, giving us new coefficients for $B$ and $C$. Thus, we obtain: $\alpha = -\frac{dy^2}{dx^2}$, $\beta = dy^2(\frac{2}{dx^2} + \frac{2}{dy^2})$, $\gamma = -1$.

The resolution using **the Fourier tridiagonal method** follows these steps:

### Domain discretization

The study domain $\Omega$ is divided into a grid of $(M+2)$ x $(N+2)$ regularly spaced points.

### Calculation of eigenvalues

We calculate the eigenvalues $d_i$ and $\delta_i$ respectively of matrices $B$ and $C$. They can be obtained as follows:

$$
\begin{cases}
d_i = 2 + 2 \left(\frac{dx}{dy}\right)^2 \left(1 - \cos{\left(\frac{\pi i}{M+1}\right)}\right) \\
\delta_i = -1 \\
i=1,...,M
\end{cases}
$$

### Calculation of the right-hand side

We define a matrix $\eta$ of dimension $N$ x $M$ as follows:

$$
\begin{cases}
\eta_{i} = Q^{t} S_i \\
i=1,...,N
\end{cases}
$$

where $Q$ is the Fourier basis change matrix, allowing the diagonalization of $B$ and $C$. It consists of the eigenvectors of $B$ and $C$ denoted $\{q_j\}_{1 \leq j \leq M}$ with

$$
\begin{cases}
q_{k,j} = \sqrt{\frac{2}{M+1}} \sin{\left(\frac{\pi j l}{M+1}\right)} \\
k=1,...,M
\end{cases}
$$

Thus, we can obtain the coefficients of $\eta$ using a discrete sine transformation as follows:

$$
\begin{cases}
\eta_{i,j} = \sqrt{\frac{2}{M+1}} \sum_{l=1}^{M}{dy^2 S_{i,l} \sin{\left(\frac{\pi j l}{M+1}\right)}} \\
i=1,...,N \quad j=1,...,M
\end{cases}
$$

### Transposition of $\eta$

Transpose the matrix $\eta$ to obtain the matrix $h$ of dimension $M$ x $N$.

### Resolution of tridiagonal systems

We define a matrix $\nu$ of dimension $M$ x $N$, where the rows are obtained by solving the following systems:

$$
T_j \nu_j = h_j \quad \quad j=1,...,M
$$

$$
\text{with} \quad T_j = 
\begin{pmatrix}
d_j & \delta_j & & & \\
\delta_j & & & & \\
& & & & \\
& & & & \delta_j \\
& & & \delta_j & d_j \\
\end{pmatrix}
$$

The matrices $T_j$ being tridiagonal, we can solve the $M$ systems using the Thomas algorithm.

### Transposition of $\nu$

Transpose the matrix $\nu$ to obtain the matrix $u$ of dimension $N$ x $M$.

### Calculation of the solution

The $Y_i$ can be obtained from $u$ as follows:

$$
\begin{cases}
Y_i = Qu_i \\
i = 1,...,N
\end{cases}
$$

The coefficients are calculated using the discrete sine transformation:

$$
\begin{cases}
Y_{i,j} = \sqrt{\frac{2}{M+1}} \sum_{l=1}^{M}{u_{i,l} \sin{\left(\frac{\pi j l}{(M + 1)}\right)}} \\
i = 1,...,N \quad j = 1,...,M
\end{cases}
$$
