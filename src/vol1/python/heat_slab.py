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
from dolfin import *
from mshr import *
import numpy as np


def createStrap(
    width_1=40.0, width_2=135.0, length_1=40.0, length_2=135.0,
    thickness=1
):

    """
    sphere = Sphere(Point(0, 0, 0), 100.0)
    cone = Cylinder(Point(0, 0, 0), Point(0, 0, -1), .35, .1)

    return cone + sphere
    """

    box_a = createBox(size_x=width_1, size_y=length_1, size_z=thickness)
    """
    box_b = createBox(
        x0=width_1, y0=-0.5*length_1,
        size_x=width_2, size_y=length_2, size_z=thickness
    )
    """

    return box_a # + box_b


def createBox(
    x0=0.0, y0=0.0, z0=0.0,
    size_x=1.0, size_y=1.0, size_z=1.0
):
    x1 = x0 + size_x
    y1 = y0 + size_y
    z1 = z0 + size_z

    print(f"({x0}, {y0}, {z0}), ({x1}, {y1}, {z1})")
    return Box(Point(x0, y0, z0), Point(x1, y1, z1))


T           = 32.0          # final time
num_steps   = 160           # number of time steps
dt          = T / num_steps # time step size
alpha       = 0.3           # parameter alpha
beta        = 0.0625        # parameter beta
theta       = 0.125         # parameter theta

msh_file = File('heat/slab_mesh.pvd')
sol_file = File('heat/slab_solution.pvd')

# Create mesh and define function space
# mesh = BoxMesh(Point(0., 0., 0.), Point(1., 0.1, 0.04), 60, 10, 5)
strap = createStrap(thickness=0.25)

meshing_domain = CSGCGALDomain3D(strap)
meshing_domain.remove_degenerate_facets(1e-12)
print(f"Meshing domain created")

gen = CSGCGALMeshGenerator3D()
gen.parameters["facet_angle"] = 25.0
gen.parameters["facet_size"] = 1.0
gen.parameters["edge_size"] = 1.0

print(f"Started meshing")
mesh = gen.generate(meshing_domain)
print(f"Finished meshing")

# mesh = generate_mesh(strap, 16, "cgal")
msh_file << mesh
print(f"Mesh saved!")

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression(
    '1 + x[0]*x[0] + alpha*x[1]*x[1] + theta*x[2]*x[2] + beta*t',
    # '1 + beta*t',
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

for n in range(num_steps):

    # Update current time
    t += dt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Plot solution
    # plot(u)
    # axi.triplot(u)
    sol_file << u

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
