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
import dolfin
from dolfin import Point, File,\
    FunctionSpace, Expression, DirichletBC, Function, TrialFunction, TestFunction, Constant,\
    interpolate, solve, errornorm, info,\
    dx, dot, grad, lhs, rhs
import meshio
import pygmsh
import numpy as np
import subprocess
import sys


def environment():
    os_version = str(subprocess.run(["uname", "-a"], capture_output=True).stdout)[2:-3]
    gmsh_version = str(subprocess.run(["gmsh", "--version"], capture_output=True).stderr)[2:-3]

    print(f"OS version (uname -a) = {os_version}")
    print(f"Python version = {sys.version}, (Debian official repository)")
    print(f"Numpy version = {np.__version__}, (Pypy repository)")
    print(f"pygmsh version = {pygmsh.__version__}, (Pypy repository: 6.1.1)")
    print(f"gmsh version = {gmsh_version}, (Debian official repository)")


def createStrap(
    width_1=40.0, length_1=40.0,
    width_2=135.0, length_2=135.0,
    thickness=1.0,
    max_size=1.0, min_size=1.0,
    part_sep=0.0001
):

    geom = pygmsh.opencascade.Geometry(
        characteristic_length_min=min_size,
        characteristic_length_max=max_size
    )

    channel_pos = [
        width_1-part_sep, -0.5*np.abs(length_1-length_2), 0.
    ]

    pad2D = geom.add_rectangle([0., 0., 0.], width_1, length_1, corner_radius=0.1)
    channel2D = geom.add_rectangle(channel_pos, width_2, length_2, corner_radius=0.1)

    strap2D = geom.boolean_union([pad2D, channel2D])
    geom.extrude(strap2D, [0., 0., thickness])

    return geom


def createMesh(object, print_info=False):

    print(f"Started meshing")
    mesh = pygmsh.generate_mesh(strap, verbose=True, extra_gmsh_arguments=None)
    mesh.prune()
    print(f"Finished meshing")
    meshio.write(meshfp, mesh)

    mesh = dolfin.Mesh()
    with dolfin.XDMFFile(meshfp) as f:
        f.read(mesh)
    if print_info:
        dolfin.info(mesh)

    return mesh


LOG_LEVEL   = 1             # 0 .. 4 (dolfin)
T           = 32.0          # final time
num_steps   = 160           # number of time steps
dt          = T / num_steps # time step size
alpha       = 0.300         # parameter alpha
beta        = 750           # parameter beta
theta       = 0.125         # parameter theta
meshfp      = '.temp/heat-slab-pygmsesh.xdmf'
solfile     = File('.temp/heat-slab-solution.pvd')

if not LOG_LEVEL or LOG_LEVEL == 0:
    dolfin.set_log_active(False)
else:
    dolfin.set_log_active(True)
dolfin.set_log_level(LOG_LEVEL)

environment()

# Create mesh and define function space
strap = createStrap(
    width_1=40.0, length_1=40.0,
    width_2=100.0, length_2=135.0,
    thickness=0.200,
    max_size=5.0, min_size=0.005
)
mesh = createMesh(strap)

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

for n in range(num_steps):

    # Update current time
    t += dt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Plot solution
    # plot(u)
    # axi.triplot(u)
    solfile << u

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
