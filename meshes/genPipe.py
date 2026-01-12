import gmsh
import math



def genPipe(r, d, l1, l2,  nr):
    """
    wenn nr nie angegeben wird, kommt es dazu, 
    dass gmsh alle weiteren Pipes mit dem selben namen erzeugt und damit
    im Endeffekt nur eine Pipe ausgibt
    """

    # meshing constraints
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.15) # 0.15
    # d > r damit revolution funktioniert

    circ1 = gmsh.model.occ.addDisk(d, 0, 0, r, r)
    circ2 = gmsh.model.occ.addDisk(-d, 0, 0, r, r)
    bend1 = gmsh.model.occ.revolve([(2, circ1)], d,d,0, 2*r, 0, 0, math.pi/2 )
    bend2 =  gmsh.model.occ.revolve([(2, circ2)], d, d,0, -2*r, 0, 0, math.pi/2)
    gmsh.model.occ.synchronize()

    v1 = [i for i in bend1 if i[0] == 3]
    v2 = [i for i in bend2 if i[0] == 3]

    gmsh.model.occ.translate(v1, -d, 0, -d)
    gmsh.model.occ.translate(v2, d, 0, d)
    print("created bends")


    gmsh.model.occ.synchronize()

    cyl1 = gmsh.model.occ.addCylinder(0, 0, -d, 0, 0, 2*d, r)
    cyl2 =  gmsh.model.occ.addCylinder(0, d, 2*d, 0, l1, 0, r)
    cyl3 =  gmsh.model.occ.addCylinder(0, d, -2*d, 0, l2, 0, r)
    cyls = [(3, c) for c in [cyl1, cyl2, cyl3]]
    gmsh.model.occ.synchronize()
    pipe = gmsh.model.occ.fuse(cyls, v1 + v2)[0]

    gmsh.model.occ.synchronize()

    #surfs = gmsh.model.getBoundary(pipe)
    #surftags = [i[1] for i in surfs if i[0] == 2]
    #gmsh.model.occ.synchronize()
    #gmsh.model.addPhysicalGroup(3, [pipe[0][1]] ,name= f"Pipe {nr}")
    #gmsh.model.addPhysicalGroup(2, surftags,name= f"Pipe {nr}")


    return pipe

if __name__ == '__main__':
    gmsh.initialize()
    print("successfully loaded gmsh lib")
    model_name = "HeatPipe"
    gmsh.model.add(model_name)
    pipe = genPipe(1, 4, 10,5,  1)
    pipe2 = genPipe(0.5, 3, 5,5,  2)
    #box = gmsh.model.occ.addBox(-1, 12 ,-9 , 2, 2, 18)
    #gmsh.model.occ.translate(pipe2, 0, 5, 0)
    gmsh.model.occ.synchronize()
    #gmsh.model.occ.fuse([(3, box)], pipe+pipe2)
    #gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.model.occ.synchronize()


    gmsh.write(f"{model_name}.msh")
    gmsh.fltk.run() # open in gmsh gui to look
    gmsh.finalize()

