import firedrake as fd
import math

###
# tags
# 1 Embedding VOl
# 5, 6, 8, 10  Outer Rand links, unten, oben, rechts

# 2 Chip TR Vol

# 3 CHIP TL VOl

# 4 Cooling Volume
# 49 Deckeplatte

###

def processFunction(t, tmax=0., sigma=1, C=10., offset=0.):
   o = C * math.exp(-1/2 * (t-tmax)**2 / sigma**2) + offset
   return fd.Constant(o)

filename = "ProtoChipCooled.pvd"
mesh = fd.Mesh("../meshes/ProtoChipCooled.msh")


print("setting up function space and functions, ... ", end="")

V = fd.FunctionSpace(mesh, "CG", 1) # testfunktionen, CG = continuous galerkin, Hutbasis

v = fd.TestFunction(V) # 
cp1 = fd.TrialFunction(V) # nach zeitschritt

# use this for c0, set up initial c distribution c(t=0, x)
VD = fd.FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten
c0 = fd.Function(V)

x, y, z = fd.SpatialCoordinate(mesh)
Cw = fd.Constant(20.) # Wassertemp
CHIPT0 = fd.Constant(40.)
alpha = fd.Constant(1) # Koeffizient der Konvektion, kann auch Ortsabh. sein

# setze Anfangstemperatur
c0.interpolate(CHIPT0, subset=mesh.cell_subset(1)) 
c0.interpolate(CHIPT0, subset=mesh.cell_subset(2)) 
c0.interpolate(CHIPT0, subset=mesh.cell_subset(3)) 
c0.interpolate(CHIPT0, subset=mesh.cell_subset(4)) 
c0.interpolate(CHIPT0, subset=mesh.cell_subset(5)) 
c0.interpolate(CHIPT0, subset=mesh.cell_subset(6)) 
c0.interpolate(CHIPT0, subset=mesh.cell_subset(7)) 


# set Diffusionskoeffizient 
D = fd.Function(VD)
D.interpolate(fd.Constant(0.3),  subset=mesh.cell_subset(1))
D.interpolate(fd.Constant(0.4),  subset=mesh.cell_subset(2))
D.interpolate(fd.Constant(0.4),  subset=mesh.cell_subset(3))
D.interpolate(fd.Constant(0.4),  subset=mesh.cell_subset(4))
D.interpolate(fd.Constant(1.),  subset=mesh.cell_subset(5))
D.interpolate(fd.Constant(1.),  subset=mesh.cell_subset(6))
D.interpolate(fd.Constant(0.5),  subset=mesh.cell_subset(7))

outfile = fd.VTKFile(filename)
print("done")


# chip BCS
bcs = [                                
        #fd.DirichletBC(V, fd.Constant(0.), [9, 10, 11])
        ]
# wir setzen adiabatische Randbedingungen, bis auf einen konvektiven Rand
t = 0.
dt = 0.05
h = fd.Constant(dt) # Zeitschritt für implizite Euler Integration
dx = fd.Measure("dx", domain=mesh)



for i in range(300):

        
    print(f"iteration: {i}", end="\r")
    # linke seite, löse für Zeitschritt
    a = cp1 * v  * dx \
            + h * D * fd.inner(fd.grad(cp1), fd.grad(v))  * dx \
            + h * D * alpha * cp1 * v  * fd.ds(8) 
    L = c0 * v  * fd.dx \
            + h * D * alpha * Cw * v  * fd.ds(8) \
            + h * v * processFunction(t+dt, tmax=13., sigma=1, C=40.) * dx(2)\
            + h * v * processFunction(t+dt, tmax=20., sigma=0.5, C=20.) * dx(3)\
            + h * v * processFunction(t+dt, tmax=7., sigma=1.2, C=30.) * dx(3)\
            + h * v * processFunction(t+dt, tmax=17., sigma=1.7, C=70.) * dx(4)


    fd.solve(a == L, c0, bcs) # speichert die Lösung cp1 in c0
    t += dt
    if i % 1 == 0:

        outfile.write(c0) 




print(f"successfully written to: {filename}")

#VV = VectorFunctionSpace(mesh, "CG", 1)
#j = project(- sigma * grad(uu), VV)
#VTKFile("j.pvd").write(j)

#VVD = VectorFunctionSpace(mesh, "DG", 0)
#jd = project(- sigma * grad(uu), VVD)
#VTKFile("jd.pvd").write(jd)

