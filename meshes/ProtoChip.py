import gmsh
import math
from genChipMesh import genChipMesh


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
model_name = "ProtoChip"
gmsh.model.add(model_name)

# meshing constraints
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.15) # 0.15

embedding = gmsh.model.occ.addBox(-2, -2, 0., 4, 4, 0.5)


chip1 = genChipMesh(1, 1, 0.2)
chip2 = gmsh.model.occ.copy(chip1)
gmsh.model.occ.synchronize()
gmsh.model.occ.translate(chip1, 1, 1, 0.5)
gmsh.model.occ.translate(chip2, -1, 1, 0.5)

cooling = genChipMesh(2, 1, 0.2)
gmsh.model.occ.translate(cooling, 0, -1, 0.5)
objs = [c[0] for c in [chip1, chip2, cooling]]
embedding = gmsh.model.occ.cut([(3, embedding)], objs, removeTool=False)[0]
vols = gmsh.model.occ.fragment(embedding, objs)[0]


gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.model.occ.synchronize()

entities = vols
surfs= gmsh.model.getBoundary(vols)
#cooling_surf = [i[1][0] for i in cooling_surf]
entities += surfs
extractObjectsIntoPhysGroups(3, entities)
extractObjectsIntoPhysGroups(2, entities)
gmsh.model.occ.synchronize()

gmsh.write(f"{model_name}.msh")
gmsh.fltk.run() # open in gmsh gui to look
gmsh.finalize()

