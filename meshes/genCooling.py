import gmsh
import math

def genCoolingBox(boxlen, boxthick, l1, dist):

    box = gmsh.model.occ.addBox(-boxlen/2, -boxlen/2, -boxthick/2, boxlen, boxlen, boxthick)

    cutout_v1 = gmsh.model.occ.addBox(-l1/2, -l1/2, -boxthick/2, l1, l1, boxthick)
    cutout_v2 = gmsh.model.occ.copy([(3, cutout_v1)])
    cutout_h = gmsh.model.occ.addBox(-boxlen*0.9/2, -l1/2, 0, boxlen*0.9, l1, l1)


    gmsh.model.occ.synchronize()

    gmsh.model.occ.translate([(3, cutout_v1)], (-boxlen+l1)/2, -boxlen/2+l1, 0)
    gmsh.model.occ.translate(cutout_v2, (boxlen-l1)/2, -boxlen/2+l1, 0)
    gmsh.model.occ.translate([(3, cutout_h)], 0, -boxlen/2+l1, boxthick/2-l1)
    bogen = gmsh.model.occ.fuse([(3, cutout_v1)], cutout_v2 + [(3, cutout_h)])[0]

    gmsh.model.occ.synchronize()

    tmp = [(3, box)]
    i = int(boxlen/(l1+dist))
    for j in range(i+1):
        if j == i:
            o = gmsh.model.occ.cut(tmp, bogen)[0]
        else:
            o = gmsh.model.occ.cut(tmp, bogen, removeTool=False)[0]
            gmsh.model.occ.translate(bogen, 0, dist+l1, 0)
            tmp = o
    return o 


if __name__ == '__main__':
    gmsh.initialize()
    print("successfully loaded gmsh lib")
    model_name = "CoolingUnit"
    gmsh.model.add(model_name)

    # meshing constraints
    #gmsh.option.setNumber("Mesh.MeshSizeMax", 0.15) # 0.15



    boxlen = 2.
    boxthick = 0.5
    l1 = 0.1
    dist = 0.2

    genCoolingBox(boxlen, boxthick, l1, dist)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.generate(3)
    gmsh.model.occ.synchronize()

    gmsh.write(f"{model_name}.msh")
    gmsh.fltk.run() # open in gmsh gui to look
    gmsh.finalize()

