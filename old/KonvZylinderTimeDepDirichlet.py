import firedrake as fd
import math

###
# tags
# 1 VOlumen
# 2 Fläche im Wasser
# 3 Fläche im Wasser
# 4 Randfläche negative z-Richtung
# 5 Fläche außerhalb des Wassers
# 6 Randfläche positive z-Richtung
###

filename = "KonvZylinderTimeDepDirichlet.pvd"
mesh = fd.Mesh("../meshes/cylinder.msh")

x, y, z = fd.SpatialCoordinate(mesh)
r = fd.sqrt(x**2 + y**2)

print("setting up function space and functions, ... ", end="")

V = fd.FunctionSpace(mesh, "CG", 1) # testfunktionen, CG = continuous galerkin, Hutbasis

v = fd.TestFunction(V) # 
cp1 = fd.TrialFunction(V) # nach zeitschritt

# use this for c0, set up initial c distribution c(t=0, x)
VD = fd.FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten
c0 = fd.Function(V)
# setze Anfangstemperatur auf 10 im gesamten Volumen
c0.interpolate(fd.Constant(0.), subset=mesh.cell_subset(1))


# set Diffusionskoeffizient = 0.1 überall
D = fd.Function(VD)
D.interpolate(fd.Constant(0.6),  subset=mesh.cell_subset(1))


uu = fd.Function(V) # output function
outfile = fd.VTKFile(filename)
outfile.write(c0) # to capture the initial conditions as well
print("done")
f = lambda t : abs(10*math.sin(t))
t = 0
C = fd.Constant(f(t))
bcs = [                                
        fd.DirichletBC(V, C, [4])
        ]

h = fd.Constant(0.2) # Zeitschritt für implizite Euler Integration
cw = fd.Constant(-5.) # Wassertemperatur
for i in range(100):
        print(f"iteration: {i}", end="\r")
        # linke seite, löse für Zeitschritt
        a = cp1 * v  * fd.dx \
                + h * D * fd.inner(fd.grad(cp1), fd.grad(v))  * fd.dx
        # hier ist schon die Randbedingung dabei, das j am Rand (3) einen konstanten Strom von 1 hat
        L = c0 * v  * fd.dx \
                - (c0 - cw) * h * v  * fd.ds(2)\
                - (c0 - cw) * h * v  * fd.ds(3) # im Wasser

        
        fd.solve(a == L, c0, bcs) # speichert die Lösung cp1 in c0
        t += 0.2
        C.assign(f(t))
        outfile.write(c0) # damit später Animationen gemacht werden können
        #c0 = uu




print(f"successfully written to: {filename}")

#VV = VectorFunctionSpace(mesh, "CG", 1)
#j = project(- sigma * grad(uu), VV)
#VTKFile("j.pvd").write(j)

#VVD = VectorFunctionSpace(mesh, "DG", 0)
#jd = project(- sigma * grad(uu), VVD)
#VTKFile("jd.pvd").write(jd)

