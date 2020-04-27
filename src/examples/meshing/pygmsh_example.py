import sys
import subprocess

import pygmsh
import meshio
import numpy as np


def setup_info():
    os_version = str(subprocess.run(["uname", "-a"], capture_output=True).stdout)[2:-3]
    gmsh_version = str(subprocess.run(["gmsh", "--version"], capture_output=True).stderr)[2:-3]

    print(f"OS version (uname -a) = {os_version}")
    print(f"Python version = {sys.version}, (Debian official repository)")
    print(f"Numpy version = {np.__version__}, (Pypy repository)")
    print(f"pygmsh version = {pygmsh.__version__}, (Pypy repository: 6.1.1)")
    print(f"gmsh version = {gmsh_version}, (Debian official repository)")


def drawCross():
    # Draw a cross.
    geom = pygmsh.built_in.Geometry()

    cross = geom.add_polygon([
            [ 0.0,  0.5, 0.0], [-0.1,  0.1, 0.0], [-0.5,  0.0, 0.0], [-0.1, -0.1, 0.0],
            [ 0.0, -0.5, 0.0], [ 0.1, -0.1, 0.0], [ 0.5,  0.0, 0.0], [ 0.1,  0.1, 0.0]
        ],
        lcar=0.05
    )

    axis = [0, 0, 1]
    angle = (2.0 / 6.0 * np.pi)

    geom.extrude(
        cross,
        translation_axis=axis,
        rotation_axis=axis,
        point_on_axis=[0, 0, 0],
        angle=angle
    )

    return geom


setup_info()

cross = drawCross()
mesh = pygmsh.generate_mesh(cross)
print(f"Geometry to mesh process finished, mesh = {mesh}")
points, cells, point_data, cell_data, field_data = mesh

meshio.write(".build/pygsm-example.xdmf", meshio.Mesh(
    points=points, cells={"triangle": cells["triangle"]})
)
