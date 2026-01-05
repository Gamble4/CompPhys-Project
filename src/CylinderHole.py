import firedrake as fd
import math

###
# tags
# 1 Outer Volume
# 6, 7, 9. 11 Outer Rand links, unten, oben, rechts
# 8, 10 Deckelplatte, Bodenplatte Outer

# 3 Embedding Volume
# 21, 22 Deckelplatte, Bodenplatte Embedding

# 2 Ring Volume
# 5, 12 Fläche Ring außen, innen
# 13, 14 Fläche Ring oben, unten

# 4 Chip Volume
# 15-18 Seitenplatten Chip
# 19, 20 Bodenplatte, Deckelplatte Chip


###


def gaussian(t, tmax, sigma, C, offset):
    o = C * math.exp(-1/2 * (t - tmax)**2 / (sigma)**2) + offset
    return fd.Constant(o)

filename = "CylinderHole.pvd"
mesh = fd.Mesh("../meshes/CylinderHole.msh")


print("setting up function space and functions, ... ", end="")

V = fd.FunctionSpace(mesh, "CG", 1) # testfunktionen, CG = continuous galerkin, Hutbasis

v = fd.TestFunction(V) # 
cp1 = fd.TrialFunction(V) # nach zeitschritt

# use this for c0, set up initial c distribution c(t=0, x)
VD = fd.FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten
c0 = fd.Function(V)

x, y, z = fd.SpatialCoordinate(mesh)
r = fd.sqrt(x**2 + y**2)
rRi = fd.Constant(1.4) # we know this from the mesh
rRo = fd.Constant(1.6)
Cw = fd.Constant(10.) # Wassertemp


PIPET0 = fd.Constant(40.)
# setze Anfangstemperatur
f = 20. * fd.cos(r/rRi * math.pi / 2) + PIPET0
g = PIPET0 - (r-rRi)/(rRo - rRi) * (PIPET0 - Cw)
c0.interpolate(Cw, subset=mesh.cell_subset(1)) 
c0.interpolate(g, subset=mesh.cell_subset(2)) 
c0.interpolate(f, subset=mesh.cell_subset(3)) 
c0.interpolate(g, subset=mesh.cell_subset(4)) 

#c0.interpolate(EMBEDT0, subset=mesh.cell_subset(5)) 

# set Diffusionskoeffizient 
D = fd.Function(VD)
D.interpolate(fd.Constant(0.4),  subset=mesh.cell_subset(1))
D.interpolate(fd.Constant(0.005),  subset=mesh.cell_subset(2))
D.interpolate(fd.Constant(0.4),  subset=mesh.cell_subset(3))
D.interpolate(fd.Constant(1.),  subset=mesh.cell_subset(4))

outfile = fd.VTKFile(filename)
outfile.write(c0) # to capture the initial conditions as well
print("done")


# chip BCS
bcs = [                                
        fd.DirichletBC(V, fd.Constant(Cw), [5, 6, 8]),
        #        fd.DirichletBC(V, f, [20, 22]),
        ]
t = 0.
dt = 0.1
h = fd.Constant(dt) # Zeitschritt für implizite Euler Integration
dx = fd.Measure("dx", domain=mesh)



for i in range(300):

        
    print(f"iteration: {i}", end="\r")
    # linke seite, löse für Zeitschritt
    a = cp1 * v  * dx \
            + h * D * fd.inner(fd.grad(cp1), fd.grad(v))  * dx \
            + D * (cp1) * h * v  * fd.ds(12) 
    L = c0 * v  * fd.dx \
            +  Cw * h * v  * fd.ds(11) \
            + h * v * f * gaussian(t+dt, 15, sigma=2., C=0.2, offset=0.) * dx(3) 

    fd.solve(a == L, c0, bcs) # speichert die Lösung cp1 in c0
    outfile.write(c0) # damit später Animationen gemacht werden können
    t += dt




print(f"successfully written to: {filename}")
