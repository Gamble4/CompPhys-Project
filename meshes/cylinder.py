import gmsh

gmsh.initialize()
print("successfully loaded gmsh lib")
model_name = "cylinder"
gmsh.model.add(model_name)

# meshing constraints
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.15) # 0.15

c1 = gmsh.model.occ.addCylinder(0,0,-5, 0,0, 10, 2) 
c2 = gmsh.model.occ.addCylinder(0,0, 5, 0,0, -10, 2,  angle=3.1415)

s = gmsh.model.occ.fuse([(3, c1)], [(3, c2)])
gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)

entities = gmsh.model.getEntities()
vols = [i for i in entities if i[0] == 3]
vol_tags = [i[1] for i in vols]
gmsh.model.addPhysicalGroup(3, vol_tags,name= "volume")

surfs = [i for i in entities if i[0] == 2]
surf_tags = [i[1] for i in surfs]
for num, surf_tag in enumerate(surf_tags):
    gmsh.model.addPhysicalGroup(2, [surf_tag], name=f"surface {num+1}")
gmsh.model.occ.synchronize()

gmsh.write(f"{model_name}.msh")
gmsh.fltk.run() # open in gmsh gui to look
gmsh.finalize()

