import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def solve_heat_equation_fem(k, T, f0, f1, g, nx, nt):
    """
    Solve the heat equation using time discretization first (implicit Euler),
    followed by FEM spatial discretization with full mass matrix

    Parameters:
    -----------
    k : float
        Diffusion coefficient
    T : float
        Final time
    f0 : function
        Left boundary condition function u(0,t)
    f1 : function
        Right boundary condition function u(1,t)
    g : function
        Initial condition function u(x,0)
    nx : int
        Number of spatial elements
    nt : int
        Number of time steps

    Returns:
    --------
    x : ndarray
        Spatial grid points
    t : ndarray
        Time grid points
    u : ndarray
        Solution u(x,t) at all grid points
    """
    # spatial and time meshes
    x = np.linspace(0, 1, nx+1)
    dx = 1.0 / nx
    t = np.linspace(0, T, nt+1)
    dt = T / nt

    # solution array, in code as u[t,x]
    u = np.zeros((nt+1, nx+1))

    # initial condition u(x,0) = g(x)
    for i in range(nx+1):
        u[0, i] = g(x[i])

    # mass matrix
    M = np.zeros((nx+1, nx+1))
    # stiffness matrix
    K = np.zeros((nx+1, nx+1))

    # element mass matrix
    m_elem = np.array([
        [2.0, 1.0],
        [1.0, 2.0]
    ]) * dx / 6.0

    # Element stiffness matrix
    k_elem = np.array([
        [1.0, -1.0],
        [-1.0, 1.0]
    ]) / dx

    # assemble global matrtix from element contributions
    for e in range(nx):
        i = e      # Left node
        j = e + 1  # Right node

        M[i, i] += m_elem[0, 0]
        M[i, j] += m_elem[0, 1]
        M[j, i] += m_elem[1, 0]
        M[j, j] += m_elem[1, 1]

        K[i, i] += k_elem[0, 0]
        K[i, j] += k_elem[0, 1]
        K[j, i] += k_elem[1, 0]
        K[j, j] += k_elem[1, 1]

    K = dt * k * K
    A = M + K

    for n in range(nt):
        # boundary conditions for the next time step
        u[n+1, 0] = f0(t[n+1])
        u[n+1, nx] = f1(t[n+1])

        # right-hand side vector: M*u^n
        b = np.zeros(nx+1)

        # compute M*u^n for interior nodes
        for i in range(1, nx):
            for j in range(nx+1):
                b[i] += M[i, j] * u[n, j]

            #
            b[i] -= K[i, 0] * u[n+1, 0]  # Left boundary contribution
            b[i] -= K[i, nx] * u[n+1, nx]  # Right boundary contribution

            # mass matrix contributions (only for nodes adjacent to boundaries)
            if i == 1:
                b[i] -= M[i, 0] * u[n+1, 0]
            if i == nx-1:
                b[i] -= M[i, nx] * u[n+1, nx]

        # solve for interior nodes
        A_interior = A[1:nx, 1:nx]
        b_interior = b[1:nx]
        u_interior = np.linalg.solve(A_interior, b_interior)

        # update u
        u[n+1, 1:nx] = u_interior

    return x, t, u

def plot_solution(x, t, u):
    """
    Plot the solution u(x,t)

    Parameters:
    -----------
    x : ndarray
        Spatial grid points
    t : ndarray
        Time grid points
    u : ndarray
        Solution u(x,t) at all grid points
    """
    X, T = np.meshgrid(x, t)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, T, u, cmap=cm.viridis, linewidth=0, antialiased=True)

    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u(x,t)')
    ax.set_title('Solution of Heat Equation using FEM')
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.savefig('heat_equation_solution.png', dpi=300)
    plt.show()

    # plot solution at different time steps
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs = axs.flatten()

    time_indices = [0, len(t)//4, len(t)//2, -1]
    for i, idx in enumerate(time_indices):
        axs[i].plot(x, u[idx, :], 'b-')
        axs[i].set_xlabel('x')
        axs[i].set_ylabel('u(x,t)')
        axs[i].set_title(f't = {t[idx]:.4f}')
        axs[i].grid(True)

    plt.tight_layout()
    plt.savefig('heat_equation_time_slices.png', dpi=300)
    plt.show()
