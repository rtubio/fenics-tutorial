import dolfin
import meshio
import pygmsh
from pprint import pprint

mesh_fp = ".build/pygmsh-example-2.xdmf"

geometry = pygmsh.built_in.Geometry()
square = geometry.add_rectangle(0, 1, 0, 1, 0)
# Deprecated method
# surface = geometry.add_physical_surface(square.surface,label=1)
geometry.add_physical(square.surface, label=1)
mesh_data = pygmsh.generate_mesh(geometry)

# For some reason, the example taken from the Internet was wrong, the
# object returned as 'mesh_data' is not a tuple, but an object with fields
#
#           points, cells, point_data, cell_data, field_data = mesh_data
#
# The previous line of code is wrong due to the stated reason.

"""
pprint(vars(mesh_data))
print(f"mesh_data.points = {mesh_data.points}")
print(f"mesh_data.cells = {mesh_data.cells}")
print(f"mesh_data.point_data = {mesh_data.point_data}")
print(f"mesh_data.cell_data = {mesh_data.cell_data}")
print(f"mesh_data.cell_data = {mesh_data.field_data}")
"""

# meshio.write(mesh_fp, meshio.Mesh(points=mesh_data.points, cells=mesh_data.cells))
meshio.write(mesh_fp, mesh_data)

mesh = dolfin.Mesh()
with dolfin.XDMFFile(mesh_fp) as infile:
    infile.read(mesh)

dolfin.info(mesh)

# Runtime ERROR is thrown with the following description:
"""
Traceback (most recent call last):
  File "vol1/python/meshing/pygms_example_2.py", line 39, in <module>
    facet.normal(0)
RuntimeError: 

*** -------------------------------------------------------------------------
*** DOLFIN encountered an error. If you are not able to resolve this issue
*** using the information listed below, you can ask for help at
***
***     fenics-support@googlegroups.com
***
*** Remember to include the error message listed below and, if possible,
*** include a *minimal* running example to reproduce the error.
***
*** -------------------------------------------------------------------------
*** Error:   Unable to find normal.
*** Reason:  Normal vector is not defined in dimension 3 (only defined when the triangle is in R^2.
*** Where:   This error was encountered inside TriangleCell.cpp.
*** Process: 0
*** 
*** DOLFIN version: 2019.1.0
*** Git changeset:  74d7efe1e84d65e9433fd96c50f1d278fa3e3f3f
*** -------------------------------------------------------------------------
"""
#
# for facet in dolfin.facets(mesh):
#     facet.normal(0)
