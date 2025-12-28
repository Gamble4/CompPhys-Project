import firedrake as fd
import numpy as np

###
# tags
# 1 Embedding Volumen
# 2 Chip Volumen
# 3-6 Seitenränder Chip Surfaces
# 7 Bodenplatte Chip Surface
# 8 Deckelplatte Chip Surface
# 9, 10, 12, 14 Seitenplatten Embedding Surface
# 11 Deckelplatte mit Loch Embedding Surface
# 13 Bodenplatte Embedding Surface
###

filename = "ProtoChip.pvd"
mesh = fd.Mesh("../meshes/ProtoChip.msh")

x, y, z = fd.SpatialCoordinate(mesh)

print("setting up function space and functions, ... ", end="")

V = fd.FunctionSpace(mesh, "CG", 1) # testfunktionen, CG = continuous galerkin, Hutbasis

v = fd.TestFunction(V) # 
cp1 = fd.TrialFunction(V) # nach zeitschritt

# use this for c0, set up initial c distribution c(t=0, x)
VD = fd.FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten
c0 = fd.Function(V)


CHIPT0 = 20.
EMBEDT0 = 10.

def processFunction(t, T0, dT, TMAX = 50., TMIN=CHIPT0):
    g = T0 + t * dT
    if g > TMAX:
        return TMAX
    elif g < TMIN:
        return TMIN
    return g

# setze Anfangstemperatur
c0.interpolate(fd.Constant(EMBEDT0), subset=mesh.cell_subset(1)) # embedding
c0.interpolate(fd.Constant(CHIPT0), subset=mesh.cell_subset(2)) # chip


# set Diffusionskoeffizient 
D = fd.Function(VD)
D.interpolate(fd.Constant(0.6),  subset=mesh.cell_subset(1))
D.interpolate(fd.Constant(0.8),  subset=mesh.cell_subset(2))

uu = fd.Function(V) # output function
outfile = fd.VTKFile(filename)
outfile.write(c0) # to capture the initial conditions as well
print("done")


C = fd.Constant(20.)
bcs = [                                
        fd.DirichletBC(V, C, [3, 4, 5, 6, 7, 8 ])
        ]
t = 0.
dt = 0.1
h = fd.Constant(dt) # Zeitschritt für implizite Euler Integration
cw = fd.Constant(5.) # Wassertemperatur

for i in range(200):

    # Temperaturanstieg in diesem Zeitfenster -> Prozess
    if t > 5 and t < 10:
        T = processFunction(i-50, CHIPT0, 0.1)
        C.assign(T)
        c0.interpolate(C, subset=mesh.cell_subset(2))
    elif t >= 10:
        C.assign(processFunction(i-100, T, -0.1))
        c0.interpolate(C, subset=mesh.cell_subset(2))
        
    print(f"iteration: {i}", end="\r")
    # linke seite, löse für Zeitschritt
    a = cp1 * v  * fd.dx \
            + h * D * fd.inner(fd.grad(cp1), fd.grad(v))  * fd.dx
    # hier ist schon die Randbedingung dabei, das j am Rand (3) einen konstanten Strom von 1 hat
    L = c0 * v  * fd.dx \
            - (c0 - cw) * h * v  * fd.ds(9)  \
            - (c0 - cw) * h * v  * fd.ds(10) 
           # - (c0 - cw) * h * v  * fd.ds(12) \
           # - (c0 - cw) * h * v  * fd.ds(14) # im Wasser

    
    fd.solve(a == L, c0, bcs) # speichert die Lösung cp1 in c0
    outfile.write(c0) # damit später Animationen gemacht werden können
    t += dt




print(f"successfully written to: {filename}")

#VV = VectorFunctionSpace(mesh, "CG", 1)
#j = project(- sigma * grad(uu), VV)
#VTKFile("j.pvd").write(j)

#VVD = VectorFunctionSpace(mesh, "DG", 0)
#jd = project(- sigma * grad(uu), VVD)
#VTKFile("jd.pvd").write(jd)

