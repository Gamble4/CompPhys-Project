import gmsh
import math
from genChipMesh import genChipMesh
from genPipe import genPipe
from genCooling import genCoolingBox 


def extractBoundarySurface(vol, name, all_in_one=False):

    vol_surf = gmsh.model.getBoundary(vol)
    if all_in_one:
        vol_surf_tags = [i[1] for i in vol_surf]
        gmsh.model.addPhysicalGroup(2, vol_surf_tags, name=name)
    else:
        i = 0
        for dim, tag in vol_surf:
            gmsh.model.addPhysicalGroup(dim, [tag], name=f"{name}{i}")
            i += 1

    return None

gmsh.initialize()
print("successfully loaded gmsh lib")
model_name = "ProtoChipCooled"
gmsh.model.add(model_name)

# meshing constraints
#gmsh.option.setNumber("Mesh.MeshSizeMax", 0.15) # 0.15

cooling_sep = 0.1

boxlen = 2.
boxthick = 0.5
chipthick = 0.2

piper = boxthick/4
d = piper + 0.01 + cooling_sep/2
l1 = 0.1
l2 = 0.1    


embedding = gmsh.model.occ.addBox(-boxlen, -boxlen, 0.,2*boxlen, 2*boxlen, boxthick)


chip1 = genChipMesh(boxlen/2, boxlen/2, chipthick)
chip2 = gmsh.model.occ.copy(chip1)
chip3 = genChipMesh(boxlen, boxlen/2, chipthick)

gmsh.model.occ.translate(chip1, boxlen/2, boxlen/2, boxthick)
gmsh.model.occ.translate(chip2, -boxlen/2, boxlen/2, boxthick)
gmsh.model.occ.translate(chip3, 0, -boxlen/2, boxthick)

gmsh.model.occ.synchronize()

pipe1 = genPipe(piper, d, l1, l2, 1)
gmsh.model.occ.rotate(pipe1, 0, 0, -boxlen, 0, 0, boxlen, angle=math.pi)
gmsh.model.occ.translate(pipe1, -boxlen/3, boxlen+l1+d,-boxthick/4)
gmsh.model.occ.synchronize()

pipe2 = gmsh.model.occ.copy(pipe1)
gmsh.model.occ.translate(pipe2, 2*boxlen/3, 0, 0)
gmsh.model.occ.synchronize()

ydiff = boxlen+l1-l2
cooling = genCoolingBox(boxlen*2, boxthick, l1, 0.3)
gmsh.model.occ.translate(cooling, 0, 0, boxthick/2)
gmsh.model.occ.rotate(cooling, -boxlen, 0, 0, 2*boxlen, 0, 0, angle=math.pi)
gmsh.model.occ.translate(cooling, 0, 0, -d)
gmsh.model.occ.synchronize()




objs = [c[0] for c in [chip1, chip2, chip3]]
embedding = gmsh.model.occ.cut([(3, embedding)], objs, removeTool=False)[0]
prev_vols = gmsh.model.getEntities(3)
print(prev_vols)
vols, vols_map = gmsh.model.occ.fragment(
        embedding,
        objs + pipe1 + pipe2 + cooling )


gmsh.model.occ.synchronize()
names = ["embedding", "chip1", "chip2", "chip3", 
        "pipe1", "pipe2", "cooling"]

ii = 0
for i, v in vols:
    gmsh.model.addPhysicalGroup(i, [v], name=names[ii])
    ii += 1
gmsh.model.occ.synchronize()


gmsh.model.occ.synchronize()

cooling_surf_tags = [i[1] for i in gmsh.model.getBoundary(cooling)]
pipe1_surf_tags = [i[1] for i in gmsh.model.getBoundary(pipe1)]
pipe2_surf_tags = [i[1] for i in gmsh.model.getBoundary(pipe2)]

unit_surf_tags  = [i[1] for i in gmsh.model.getBoundary(embedding)]
chip1_surf_tags = [i[1] for i in gmsh.model.getBoundary(chip1)]
chip2_surf_tags = [i[1] for i in gmsh.model.getBoundary(chip2)]
chip3_surf_tags = [i[1] for i in gmsh.model.getBoundary(chip3)]

tmp = []

for tag in cooling_surf_tags:
    if tag not in pipe1_surf_tags:
        if tag not in pipe2_surf_tags:
            tmp.append(tag)


gmsh.model.addPhysicalGroup(2, tmp, name="cooling surf") 



tmp = []

for tag in chip1_surf_tags:
    if tag not in unit_surf_tags:
        tmp.append(tag)


gmsh.model.addPhysicalGroup(2, tmp, name="chip1 surf") 


tmp = []

for tag in chip2_surf_tags:
    if tag not in unit_surf_tags:
        tmp.append(tag)


gmsh.model.addPhysicalGroup(2, tmp, name="chip2 surf") 

tmp = []

for tag in chip3_surf_tags:
    if tag not in unit_surf_tags:
        tmp.append(tag)


gmsh.model.addPhysicalGroup(2, tmp, name="chip3 surf") 



gmsh.model.addPhysicalGroup(2, pipe1_surf_tags, name="pipe1 surf") 
gmsh.model.addPhysicalGroup(2, pipe2_surf_tags, name="pipe2 surf") 



gmsh.model.occ.synchronize()


gmsh.model.mesh.generate(3)
gmsh.model.occ.synchronize()

gmsh.model.occ.synchronize()

gmsh.write(f"{model_name}.msh")
gmsh.fltk.run() # open in gmsh gui to look
gmsh.finalize()

