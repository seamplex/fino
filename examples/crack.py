import gmsh
import sys

gmsh.initialize(sys.argv)

gmsh.model.add("square with cracks")

surf1 = 1
gmsh.model.occ.addRectangle(0, 0, 0, 1, 1, surf1)

pt1 = gmsh.model.occ.addPoint(0.0, 0.5, 0)
pt2 = gmsh.model.occ.addPoint(0.4, 0.5, 0)
line1 = gmsh.model.occ.addLine(pt1, pt2)

o, m = gmsh.model.occ.fragment([(2, surf1)], [(1, line1)])
gmsh.model.occ.synchronize()

# m contains, for each input entity (surf1, line1 and line2), the child entities
# (if any) after the fragmentation, as lists of tuples. To apply the crack
# plugin we group all the intersecting lines in a physical group

new_surf = m[0][0][1]
new_lines = [item[1] for sublist in m[1:] for item in sublist]

gmsh.model.addPhysicalGroup(2, [new_surf], 100)
gmsh.model.setPhysicalName (2, 100, "bulk")

gmsh.model.addPhysicalGroup(1, new_lines, 101)

gmsh.model.addPhysicalGroup(1, [7,9], 102)
gmsh.model.setPhysicalName (1, 102, "left")

gmsh.model.addPhysicalGroup(1, [8], 103)
gmsh.model.setPhysicalName (1, 103, "right")

gmsh.model.addPhysicalGroup(1, [6], 104)
gmsh.model.setPhysicalName (1, 104, "bottom")

gmsh.model.addPhysicalGroup(1, [10], 105)
gmsh.model.setPhysicalName (1, 105, "top")

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.03)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01)

gmsh.model.mesh.generate(2)

gmsh.plugin.setNumber("Crack", "PhysicalGroup", 101)
gmsh.plugin.run("Crack")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write("crack.msh")

gmsh.finalize()
