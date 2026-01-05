import gmsh


def genTriangleSurf(p1, p2, p3):
    """
    order :
        p1 : bottom
        p2 : top
        p3 : front
    """

    po1 = gmsh.model.occ.addPoint(*p1)
    po2 =  gmsh.model.occ.addPoint(*p2)
    po3 =  gmsh.model.occ.addPoint(*p3)
    l1 =  gmsh.model.occ.addLine(po1, po2)
    l2 =  gmsh.model.occ.addLine(po2, po3)
    l3 =  gmsh.model.occ.addLine(po3, po1)
    c1 = gmsh.model.occ.addCurveLoop([l1, l2, l3])
    s1 = gmsh.model.occ.addPlaneSurface([c1])
    gmsh.model.occ.synchronize()
    return s1

def genChipMesh(x, y, z):
    """
    creates the mesh around the center (0, 0, 0)
    and symmetricially around the z-plane
    params:
        x : x extent
        y : y extent 
        z : z extent
    """
    #embedding = gmsh.model.occ.addBox(-2, -2, -0.5, 4, 4, 1)
    xh = x/2
    yh = y/2
    zh = z/2

    b1 = gmsh.model.occ.addBox(-xh, -yh, -zh, x, y, z)

    p1 = (-xh, -yh, -zh)
    p2 = (-xh, -yh, zh)
    p3 = (-(xh+zh), -yh, 0.)
    s1 = genTriangleSurf(p1, p2, p3)
    gmsh.model.occ.synchronize()

    v1= gmsh.model.occ.extrude([(2, s1)], 0, y, 0)
    gmsh.model.occ.synchronize()


    p1 = (xh, -yh, -zh)
    p2 = (xh, -yh, zh)
    p3 = ((xh+zh), -yh, 0.)
    s2 = genTriangleSurf(p1, p2, p3)
    gmsh.model.occ.synchronize()
    v2 = gmsh.model.occ.extrude([(2, s2)], 0, y, 0)
    gmsh.model.occ.synchronize()


    p1 = (-xh, -yh, -zh)
    p2 = (-xh, -yh, zh)
    p3 = (-xh, -(yh+zh), 0.)
    s3 = genTriangleSurf(p1, p2, p3)
    gmsh.model.occ.synchronize()
    v3 = gmsh.model.occ.extrude([(2, s3)], x, 0, 0)
    gmsh.model.occ.synchronize()


    p1 = (-xh, yh, -zh)
    p2 = (-xh, yh, zh)
    p3 = (-xh, yh+zh, 0.)
    s4 = genTriangleSurf(p1, p2, p3)
    gmsh.model.occ.synchronize()
    v4 = gmsh.model.occ.extrude([(2, s4)], x, 0, 0)
    gmsh.model.occ.synchronize()

    vols = [i for i in v1+v2+v3+v4 if i[0] == 3]

    chip1 = gmsh.model.occ.fuse([(3, b1)], vols)[0]
    return chip1

