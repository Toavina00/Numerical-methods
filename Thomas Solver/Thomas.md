### Thomas Algorithm (Tridiagonal Matrix Algorithm)

The **Thomas algorithm**, also known as the **Tridiagonal Matrix Algorithm (TDMA)**, is an efficient method specifically designed for solving systems of linear equations where the coefficient matrix is **tridiagonal**. A tridiagonal matrix is an $N \times N$ matrix with non-zero elements only on the main diagonal, the first superdiagonal (the diagonal immediately above the main diagonal), and the first subdiagonal (the diagonal immediately below the main diagonal).

#### Steps of the Thomas Algorithm

Consider a system of $N$ linear equations represented in matrix form as $A \mathbf{x} = \mathbf{d}$, where $A$ is an $N \times N$ tridiagonal matrix, $\mathbf{x}$ is the $N \times 1$ vector of unknowns, and $\mathbf{d}$ is the $N \times 1$ vector of constants. The matrix $A$ has the following structure:

$$
A = \begin{pmatrix}
b_1 & c_1 & 0 & \cdots & 0 \\
a_2 & b_2 & c_2 & \cdots & 0 \\
0 & a_3 & b_3 & \ddots & \vdots \\
\vdots & \ddots & \ddots & \ddots & c_{N-1} \\
0 & \cdots & 0 & a_N & b_N
\end{pmatrix}
$$

Here:

* $a_i$ are the elements of the **subdiagonal**, where $i = 2, 3, \ldots, N$.
* $b_i$ are the elements of the **main diagonal**, where $i = 1, 2, \ldots, N$.
* $c_i$ are the elements of the **superdiagonal**, where $i = 1, 2, \ldots, N-1$.

The Thomas algorithm consists of two sequential phases: **forward elimination** and **backward substitution**.

**1. Forward Elimination (Forward Sweep):**

This step transforms the original tridiagonal system into an upper triangular system. It involves modifying the coefficients of the superdiagonal ($c_i$) and the right-hand side vector ($d_i$) row by row.

For $i = 1$:
$$
c'_1 = \frac{c_1}{b_1}
$$
$$
d'_1 = \frac{d_1}{b_1}
$$

For $i = 2, 3, \ldots, N$:
$$
c'_i = \begin{cases}
\frac{c_i}{b_i - a_i c'_{i-1}} & \text{for } i = 2, 3, \ldots, N-1 \\
0 & \text{for } i = N
\end{cases}
$$
$$
d'_i = \frac{d_i - a_i d'_{i-1}}{b_i - a_i c'_{i-1}}
$$

Here, $a_i$ are the subdiagonal elements, $b_i$ are the main diagonal elements, $c_i$ are the superdiagonal elements, and $d_i$ are the elements of the right-hand side vector. The primed coefficients $c'_i$ and $d'_i$ are the modified values. Note that for the last row ($i=N$), $c'_N$ is set to 0 because there is no $c_N$ in the original matrix.

**2. Backward Substitution (Back Sweep):**

Once the forward elimination is complete, the system is in an upper triangular form, which can be easily solved for the unknowns $x_i$ starting from the last equation and working backward.

For the last unknown ($i = N$):
$$
x_N = d'_N
$$

For the preceding unknowns ($i = N-1, N-2, \ldots, 1$):
$$
x_i = d'_i - c'_i x_{i+1}
$$

In essence, the Thomas algorithm is a specialized and efficient form of Gaussian elimination tailored to the structure of tridiagonal matrices. It significantly reduces the number of operations required compared to standard Gaussian elimination, making it a valuable tool for solving many scientific and engineering problems.
