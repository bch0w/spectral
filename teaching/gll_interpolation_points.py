"""
Modified from ChatGPT code prompt 
Can you write me a short python code that plots GLL node interpolation points
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval, legder

def gll_nodes(N):
    """Compute N+1 Gauss–Lobatto–Legendre (GLL) nodes in [-1, 1]."""
    if N == 1:
        return np.array([-1.0, 1.0])

    # Derivative of Legendre polynomial of degree N
    Pn = np.zeros(N+1)
    Pn[-1] = 1
    dPn = legder(Pn)

    # Initial guess (Chebyshev)
    x = -np.cos(np.pi * np.arange(N+1) / N)

    # Newton-Raphson refinement
    for _ in range(100):
        Px = legval(x, Pn)
        dPx = legval(x, dPn)
        dx = -((1 - x**2) * dPx + x * Px) / (2 * dPx**2 - Px * legval(x, legder(dPn)))
        x += dx
        if np.max(np.abs(dx)) < 1e-14:
            break

    x[0], x[-1] = -1, 1
    return np.sort(x)

def lagrange_basis(x, nodes, i):
    """Evaluate the i-th Lagrange basis polynomial at points x."""
    xi = nodes[i]
    L = np.ones_like(x)
    for j, xj in enumerate(nodes):
        if j != i:
            L *= (x - xj) / (xi - xj)
    return L

# Parameters
N = 6
nodes = gll_nodes(N)

# Evaluate on fine grid
x_plot = np.linspace(-1, 1, 400)
plt.figure(figsize=(20, 3), dpi=200)

# Plot each basis function
for i in range(N+1):
    Li = lagrange_basis(x_plot, nodes, i)
    plt.plot(x_plot, Li, label=f'$L_{i}(x)$', c="k", zorder=9)
    plt.plot(nodes[i], 0, 'ko', markersize=10, zorder=10)  # node marker

# plt.title(f'Lagrange Basis Functions for GLL Nodes (N={N})')
plt.xlabel('x')
plt.ylabel('$L_i(x)$')
plt.grid(True, linestyle='--', alpha=0.5)
# plt.legend(loc='upper right', fontsize=9)
plt.tight_layout()
plt.axis("off")
plt.savefig("gll_basis.png", transparent=True)
plt.show()

