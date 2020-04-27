# Copyright (C) 2014 Benjamin Kehlet
#
# This file is part of mshr.
#
# mshr is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mshr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with mshr. If not, see <http://www.gnu.org/licenses/>.
#

import dolfin
import pygmsh
# from mshr import Sphere, Cylinder, CSGCGALDomain3D, generate_mesh

# dolfin.set_log_level(dolfin.TRACE)
dolfin.set_log_active(True)
dolfin.set_log_level(4)

# Define 3D geometry
# sphere = Sphere(dolfin.Point(0, 0, 0), 0.5)
# cone = Cylinder(dolfin.Point(0, 0, 0), dolfin.Point(0, 0, -1), .35, .1)
# geometry = cone + sphere

geom = pygmsh.opencascade.Geometry(
    characteristic_length_min=0.1, characteristic_length_max=0.1
)

sphere_a = geom.add_ball([0., 0., 0.], 0.5)
cone = geom.add_cylinder([0., 0., 0.], [0., 0., 1.], .1)
figure = geom.boolean_union([sphere_a, cone])
mesh = pygmsh.generate_mesh(geom, verbose=True)

# Geometry surfaces can be saved to off files
# which can be viewed by eg. MeshLab
# meshing_domain = CSGCGALDomain3D(geometry)
# meshing_domain.remove_degenerate_facets(1e-12)
# meshing_domain.save("./.build/icecream-pygmsh-domain.off")

# Test printing
# dolfin.info("\nCompact output of 3D geometry:")
# dolfin.info(mesh)
# dolfin.info("\nVerbose output of 3D geometry:")
# dolfin.info(mesh, True)

# Generate and plot mesh
# m = generate_mesh(geometry, 16, "cgal")

# dolfin.info(m)
# dolfin.plot(m, "3D mesh")
# vtkfile = dolfin.File('.build/icecream-pygmsh-mesh.pvd')
# vtkfile << mesh

# dolfin.interactive()
