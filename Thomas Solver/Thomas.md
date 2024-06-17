### Thomas Algorithm (Tridiagonal Matrix Algorithm)

The **Thomas algorithm** is a specialized method for solving systems of linear equations where the coefficient matrix is tridiagonal (having nonzero elements only on the main diagonal, one diagonal above, and one diagonal below).

#### Steps of the Thomas Algorithm

1. **Decomposition**:
   - For a tridiagonal system $A \mathbf{x} = \mathbf{d}$, where $A$ is tridiagonal, decompose the matrix $A$ into:
     $$
     A = \begin{pmatrix}
     b_1 & c_1 & 0 & \cdots & 0 \\
     a_2 & b_2 & c_2 & \cdots & 0 \\
     0 & a_3 & b_3 & \cdots & 0 \\
     \vdots & \vdots & \vdots & \ddots & \vdots \\
     0 & \cdots & 0 & a_{N-1} & b_{N-1}
     \end{pmatrix}
     $$
     where $a_i, b_i, c_i$ are the elements of the subdiagonal, main diagonal, and superdiagonal respectively.

2. **Forward Sweep**:
   - Perform forward elimination to transform $A$ into an upper triangular matrix:
     $$
     c'_i = \begin{cases}
     \frac{c_i}{b_i} & \text{if } i = 1 \\
     \frac{c_i}{b_i - a_i \cdot c'_{i-1}} & \text{if } i = 2, 3, \ldots, N-1
     \end{cases}
     $$
     $$
     d'_i = \begin{cases}
     \frac{d_i}{b_i} & \text{if } i = 1 \\
     \frac{d_i - a_i \cdot d'_{i-1}}{b_i - a_i \cdot c'_{i-1}} & \text{if } i = 2, 3, \ldots, N
     \end{cases}
     $$

3. **Back Substitution**:
   - Perform back substitution to solve for $\mathbf{x}$ in $A \mathbf{x} = \mathbf{d}$:
     $$
     x_N = d'_N
     $$
     $$
     x_i = d'_i - c'_i \cdot x_{i+1} \quad \text{for } i = N-1, N-2, \ldots, 1
     $$
