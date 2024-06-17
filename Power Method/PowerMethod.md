### Power Algorithm

The **power iteration algorithm** is a numerical method used to find the dominant eigenvalue and its associated eigenvector of a square matrix $A$. It iteratively improves an initial guess for the eigenvector until convergence.

#### Steps of the Power Algorithm

1. **Initialization**:
   - Start with an initial guess for the eigenvector $x^{(0)}$. Typically, $x^{(0)}$ is chosen randomly or based on prior knowledge.

2. **Iteration**:
   - Update the vector iteratively using:
     $$
     x^{(k+1)} = A x^{(k)}
     $$
     where $x^{(k)}$ is the vector at iteration $k$ and $A$ is the matrix.

3. **Normalization**:
   - Normalize $x^{(k+1)}$ after each iteration to prevent it from growing indefinitely:
     $$
     x^{(k+1)} = \frac{x^{(k+1)}}{\| x^{(k+1)} \|_2}
     $$
     where $\| \cdot \|_2$ denotes the Euclidean norm.

4. **Convergence Check**:
   - Monitor the convergence by checking if the change in the estimated eigenvalue or the vector $x^{(k+1)}$ falls below a specified tolerance.

5. **Eigenvalue Estimation**:
   - Estimate the dominant eigenvalue $\lambda_1$ using the Rayleigh quotient:
     $$
     \lambda_1^{(k)} = \frac{(x^{(k)})^T A x^{(k)}}{(x^{(k)})^T x^{(k)}}
     $$
     This quotient converges to $\lambda_1$ as $k$ increases.