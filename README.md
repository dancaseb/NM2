# Numerické metody, projekt č. 4

This project implements a solver for the heat equation using the Finite Element Method (FEM) with time discretization using Euler's method.

## Problem Statement
Solve

$$\frac{\partial u}{\partial t} - k \frac{\partial^2 u}{\partial x^2} = 0$$

 for $x ∈ (0,1), t ∈ (0,T)$

With boundary conditions:
- $u(0,t) = f₀(t)$
- $u(1,t) = f₁(t)$

And initial condition:
- $u(x,0) = g(x)$

Where:
- $k \in \mathbb{R}^+$ (diffusion coefficient)
- $T \in \mathbb{R}^+$ (final time)
- $f₀(t) : (0,T) \rightarrow \mathbb{R}$ and f₁(t) : $(0,T) \rightarrow \mathbb{R}$ are boundary functions
- $g(x) : (0,1) \rightarrow \mathbb{R}$  is the initial condition function

## Implementation

The implementation consists of the following files:

1. `heat_equation_fem.py` - Basic implementation of the FEM solver
2. `examples.ipynb` - interactive notebook with a few examples

## Mathematical Approach to solve the heat equation

### Time Discretization
Euler's method is used for time discretization. 

### Spatial Discretization (FEM)

The Finite Element Method is used for spatial discretization. We use linear elements on a uniform mesh of the domain [0,1].

For each element, we construct local stiffness and mass matrices, which are then assembled into global matrices.



<!-- ## Examples

The `heat_equation_examples.py` file contains several examples:

1. **Example 1**: Heat equation with zero boundary conditions and sine initial condition
   - This has an analytical solution for comparison

2. **Example 2**: Heat equation with non-zero boundary conditions

3. **Example 3**: Heat equation with time-dependent boundary conditions

4. **Convergence Study**: Analysis of the error as the mesh is refined -->

<!-- ## Usage

To run the basic solver:
```python
python heat_equation_fem.py
```

To run the advanced solver with examples:
```python
python heat_equation_examples.py
``` -->

<!-- ## Requirements

- NumPy
- SciPy
- Matplotlib -->

<!-- ## Results

The solver produces visualizations of the solution u(x,t) over the domain Ω = (0,1)×(0,T), including:
- 3D surface plots
- 2D slices at different time points
- Comparison with analytical solutions where available
- Convergence analysis -->

## Steps to solve the heat equation
### 1. Time Discretization (Implicit Euler Method)

The first step is to discretize the heat equation in time using the Implicit Euler method.

### Time Grid

We divide the time interval $[0,T]$ into $n_t$ equal subintervals:
- Time step: $\Delta t = \frac{T}{n_t}$
- Time points: $t_n = n \cdot \Delta t$ for $n = 0, 1, 2, \ldots, n_t$

### Implicit Euler Discretization

Using the implicit euler method, we can approximate the time derivation:

$$\frac{\partial u}{\partial t}(x, t_{n+1}) \approx \frac{u(x, t_{n+1}) - u(x, t_n)}{\Delta t}$$

We can plug in this term into our original equation and get

$$\frac{u(x, t_{n+1}) - u(x, t_n)}{\Delta t} = k \frac{\partial^2 u}{\partial x^2}(x, t_{n+1})$$

After rearranging to isolate the terms at time $t_{n+1}$:

$$u(x, t_{n+1}) - \Delta t \cdot k \frac{\partial^2 u}{\partial x^2}(x, t_{n+1}) = u(x, t_n)$$

This gives us a sequence of problems to solve at each time step. For each $n = 0, 1, 2, \ldots, n_t-1$, we need to solve:

$$u^{n+1}(x) - \Delta t \cdot k \frac{d^2 u^{n+1}}{dx^2}(x) = u^n(x)$$

where $u^n(x) = u(x, t_n)$.

----------------------------

### 2. Weak Formulation

In the next step we will convert the strong form of the equation to its weak form (so that we can apply the FEM)

### Multiply by a Test Function

Multiply both sides of the time-discretized equation by a test function $v(x)$ and integrate over the domain $\Omega = [0,1]$:

$$\int_0^1 u^{n+1}(x) v(x) dx - \Delta t \cdot k \int_0^1 \frac{d^2 u^{n+1}}{dx^2}(x) v(x) dx = \int_0^1 u^n(x) v(x) dx$$

### Integration by Parts

Using integration by parts on the term with second order derivative we get

$$\int_0^1 \frac{d^2 u^{n+1}}{dx^2}(x) v(x) dx = \left[ \frac{du^{n+1}}{dx}(x) v(x) \right]_0^1 - \int_0^1 \frac{du^{n+1}}{dx}(x) \frac{dv}{dx}(x) dx$$


### Weak Form

The weak form becomes:

$$\int_0^1 u^{n+1}(x) v(x) dx + \Delta t \cdot k \int_0^1 \frac{du^{n+1}}{dx}(x) \frac{dv}{dx}(x) dx = \int_0^1 u^n(x) v(x) dx + \text{boundary terms}$$


----------------------------

### 3. FEM Spatial Discretization

Now we apply the FEM to discretize the weak form in space.

### Spatial Grid

We divide the spatial domain $[0,1]$ into $n_x$ equal elements:
- Element size: $\Delta x = \frac{1}{n_x}$
- Grid points: $x_i = i \cdot \Delta x$ for $i = 0, 1, 2, \ldots, n_x$

### Basis Functions

We use piecewise linear basis functions (hat functions) $\phi_i(x)$ for $i = 0, 1, 2, \ldots, n_x$:

$$\phi_i(x) = 
\begin{cases}
\frac{x - x_{i-1}}{x_i - x_{i-1}} & \text{if } x \in [x_{i-1}, x_i] \\
\frac{x_{i+1} - x}{x_{i+1} - x_i} & \text{if } x \in [x_i, x_{i+1}] \\
0 & \text{otherwise}
\end{cases}$$

### Approximation of the Solution

We approximate the solution at time $t_{n+1}$ as:

$$u^{n+1}(x) \approx \sum_{j=0}^{n_x} u^{n+1}_j \phi_j(x)$$

where $u^{n+1}_j$ are the unknown coefficients representing the solution values at the grid points.

### Galerkin Method
In the galerkin method, we choose the test functions to be the same as the basis functions $v(x) = \phi_i(x)$ for $i = 0, 1, 2, \ldots, n_x$.

Substituting into the weak form, we get:

$$\sum_{j=0}^{n_x} u^{n+1}_j \int_0^1 \phi_j(x) \phi_i(x) dx + \Delta t \cdot k \sum_{j=0}^{n_x} u^{n+1}_j \int_0^1 \frac{d\phi_j}{dx}(x) \frac{d\phi_i}{dx}(x) dx = \sum_{j=0}^{n_x} u^{n}_j \int_0^1 \phi_j(x) \phi_i(x) dx + \text{boundary terms}$$

### Matrix Form

In order to write our program, we can write the system in matrix form:

$$(M + \Delta t \cdot k \cdot K) U^{n+1} = M U^n + \text{boundary terms}$$

where:
- $M$ is the mass matrix with entries $M_{ij} = \int_0^1 \phi_j(x) \phi_i(x) dx$
- $K$ is the stiffness matrix with entries $K_{ij} = \int_0^1 \frac{d\phi_j}{dx}(x) \frac{d\phi_i}{dx}(x) dx$
- $U^n$ is the vector of solution values at time $t_n$

### Mass and Stiffness Matrices construction
For an element of length $dx$, the local mass matrix is: $$m_{elem} = \frac{dx}{6} \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}$$
For an element of length $dx$, the local stiffness matrix is: $$k_{elem} = \frac{1}{h} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$
These local matrices are then assembled into the global matrices based on the connectivity of the elements.

----------------------------

### 4. Handling Boundary Conditions

For Dirichlet boundary conditions, we set:
- $u^{n+1}_0 = f_0(t_{n+1})$
- $u^{n+1}_{n_x} = f_1(t_{n+1})$

We then solve the system only for the interior nodes, adjusting the right-hand side to account for the boundary values.

