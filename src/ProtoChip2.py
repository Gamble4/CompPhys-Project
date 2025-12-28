import firedrake as fd
import math
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

def UFLabs(f):
    return fd.conditional(f<0, -f, f)


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



TCHIP0 = 20.
TEMBED0 = 10.
# setze Anfangstemperatur
c0.interpolate(fd.Constant(TEMBED0), subset=mesh.cell_subset(1)) # embedding
c0.interpolate(fd.Constant(TCHIP0), subset=mesh.cell_subset(2)) # chip


# set Diffusionskoeffizient 
D = fd.Function(VD)
D.interpolate(fd.Constant(0.6),  subset=mesh.cell_subset(1))
D.interpolate(fd.Constant(0.8),  subset=mesh.cell_subset(2))

uu = fd.Function(V) # output function
outfile = fd.VTKFile(filename)
outfile.write(c0) # to capture the initial conditions as well
print("done")


C1 = fd.Constant(TCHIP0)
# simuliere ein Gitter, auf dem das Embedding liegt
# niedrige Temperatur -> Metall
# hohe temperatur -> Luft
L = 2.; nu = 5.; amp = 10.
f =  amp*fd.sin(x * nu *  2*np.pi/L) + amp*fd.sin(y*nu*2*np.pi/L) + 10.
C2 = UFLabs(f)


processFunction = lambda t : TCHIP0 + t * 1


bcs = [                                
        fd.DirichletBC(V, C1, [3, 4, 5, 6, 7, 8 ]),
        fd.DirichletBC(V, C2, [13])
        ]
t = 0.
dt = 0.1
h = fd.Constant(dt) # Zeitschritt für implizite Euler Integration
cw = fd.Constant(5.) # Wassertemperatur
for i in range(200):

    # Temperaturanstieg in diesem Zeitfenster -> Prozess
    if t > 5 and t < 10:
        C1.assign(40.)
        c0.interpolate(C1, subset=mesh.cell_subset(2))
    else: # Prozess vorbei -> sofortige Abkühlung, nicht realistisch
        C1.assign(TCHIP0)
        c0.interpolate(C1, subset=mesh.cell_subset(2))
        
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

