from dolfin import *

"""
mesh = UnitIntervalMesh(10)
plot(mesh, title="Unit interval")

mesh = UnitSquareMesh(10, 10)
plot(mesh, title="Unit square")

mesh = UnitSquareMesh(10, 10, "left")
plot(mesh, title="Unit square (left)")

mesh = UnitSquareMesh(10, 10, "crossed")
plot(mesh, title="Unit square (crossed)")

mesh = UnitSquareMesh(10, 10, "right/left")
plot(mesh, title="Unit square (right/left)")

mesh = RectangleMesh(0.0, 0.0, 10.0, 4.0, 10, 10)
plot(mesh, title="Rectangle")

mesh = RectangleMesh(-3.0, 2.0, 7.0, 6.0, 10, 10, "right/left")
plot(mesh, title="Rectangle (right/left)")

if has_cgal():
    mesh = CircleMesh(Point(0.0, 0.0), 1.0, 0.2)
    plot(mesh, title="Circle (unstructured)")

    mesh = EllipseMesh(Point(0.0, 0.0), [3.0, 1.0], 0.2)
    plot(mesh, title="Ellipse mesh (unstructured)")

    mesh = SphereMesh(Point(0.0, 0.0, 0.0), 1.0, 0.2)
    plot(mesh, title="Sphere mesh (unstructured)")

    mesh = EllipsoidMesh(Point(0.0, 0.0, 0.0), [3.0, 1.0, 2.0], 0.2)
    plot(mesh, title="Ellipsoid mesh (unstructured)")
"""

mesh = UnitCubeMesh(10, 10, 10)
# plot(mesh, title="Unit cube")


"""
mesh = BoxMesh(0.0, 0.0, 0.0, 10.0, 4.0, 2.0, 10, 10, 10)
plot(mesh, title="Box")
"""

interactive()