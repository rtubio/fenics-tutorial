"""
FEniCS tutorial demo program: Heat equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.

  u'= Laplace(u) + f  in the unit square
  u = u_D             on the boundary
  u = u_0             at t = 0

  u = 1 + x^2 + alpha*y^2 + beta*t
  f = beta - 2 - 2*alpha
"""

from __future__ import print_function
from dolfin import Point, File,\
    FunctionSpace, Expression, DirichletBC, Function, TrialFunction, TestFunction, Constant,\
    interpolate, solve, errornorm, info,\
    UnitCubeMesh,\
    dx, dot, grad, lhs, rhs
import numpy as np


T           = 32.0          # final time
num_steps   = 160           # number of time steps
dt          = T / num_steps # time step size
alpha       = 0.3           # parameter alpha
beta        = 0.0625        # parameter beta
theta       = 0.125         # parameter theta

# Create mesh and define function space
nx = ny = nz = 8
mesh = UnitCubeMesh(nx, ny, nz)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression(
    '1 + x[0]*x[0] + alpha*x[1]*x[1] + theta*x[2]*x[2] + beta*t',
    degree=2, alpha=alpha, theta=theta, beta=beta, t=0
)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define initial value
u_n = interpolate(u_D, V)
# u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time-stepping
u = Function(V)
t = 0
vtkfile = File('heat/solution_3D.pvd')

for n in range(num_steps):

    # Update current time
    t += dt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Plot solution
    # plot(u)
    # axi.triplot(u)
    vtkfile << u

    # Compute error at vertices
    # u_e = interpolate(u_D, V)
    # error = np.abs(u_e.vector().array() - u.vector().array()).max()
    # print('t = %.2f: error = %.3g' % (t, error))

    # Update previous solution
    u_n.assign(u)

# Compute error in L2 norm
error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Print errors
print('error_L2  =', error_L2)
print('error_max =', error_max)

# Hold plot
# interactive()
# plt.show()
