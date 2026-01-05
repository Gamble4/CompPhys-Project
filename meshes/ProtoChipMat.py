import gmsh


def extractObjectsIntoPhysGroups(dim : int, entity_list : list):
    """
    packs all entities of given dim into a unique physical group
    """
    
    match dim:
        case 0:
            s = "Point"
        case 1:
            s = "Line"
        case 2:
            s = "Surface"
        case 3:
            s = "Volume"
        case _:
            raise ValueError("only dims of 0-3 valid")

    e = [i for i in entity_list if i[0] == dim]
    e_tags = [i[1] for i in e]
    for num, tag in enumerate(e_tags):
        gmsh.model.addPhysicalGroup(dim, [tag],name= f"{s} {num+1}")
    
    return None

gmsh.initialize()
print("successfully loaded gmsh lib")
model_name = "ProtoChipMat"
gmsh.model.add(model_name)

# meshing constraints
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.15) # 0.15

b1 = gmsh.model.occ.addBox(-2, -2, -0.5, 4, 4, 1)
c1 = gmsh.model.occ.addCylinder(0, 0, -0.5, 0, 0, 1, r=1.6)
outer = gmsh.model.occ.cut([(3, b1)], [(3, c1)], removeTool=False)


c2 = gmsh.model.occ.addCylinder(0, 0, -0.5, 0, 0, 1, r=1.3)
ring = gmsh.model.occ.cut([(3, c1)], [(3, c2)], removeTool=False)[0]

chip = gmsh.model.occ.addBox(-0.5, -0.5, 0, 1, 1, 0.5)
embedding = gmsh.model.occ.cut([(3, c2)], [(3, chip)], removeTool=False)[0]

gmsh.model.occ.synchronize()


gmsh.model.mesh.generate(3)
gmsh.model.occ.synchronize()

entities = gmsh.model.getEntities()
extractObjectsIntoPhysGroups(3, entities)
extractObjectsIntoPhysGroups(2, entities)
gmsh.model.occ.synchronize()

gmsh.write(f"{model_name}.msh")
gmsh.fltk.run() # open in gmsh gui to look
gmsh.finalize()

