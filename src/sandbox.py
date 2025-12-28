import firedrake as fd

filename = "sandbox.pvd"
mesh = fd.Mesh("../meshes/bone_with_impurity.msh")

print("setting up function space and functions, ... ", end="")

V = fd.FunctionSpace(mesh, "CG", 1) # testfunktionen, CG = continuous galerkin, Hutbasis

v = fd.TestFunction(V) # 
cp1 = fd.TrialFunction(V) # nach zeitschritt

# use this for c0, set up initial c distribution c(t=0, x)
VD = fd.FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten
c0 = fd.Function(V)
c0.interpolate(fd.Constant(0.), subset=mesh.cell_subset(1))
c0.interpolate(fd.Constant(0.), subset=mesh.cell_subset(2))

# set Diffusionskoeffizient = 1 überall
D = fd.Function(VD)
D.interpolate(fd.Constant(1.),  subset=mesh.cell_subset(1))
D.interpolate(fd.Constant(1.),  subset=mesh.cell_subset(2))


uu = fd.Function(V) # output function
outfile = fd.VTKFile(filename)
outfile.write(c0) # to capture the initial conditions as well
print("done")


h = fd.Constant(0.05) # Zeitschritt für implizite Euler Integration
for i in range(20):
        print(f"iteration: {i}", end="\r")
        # linke seite, löse für Zeitschritt
        a = cp1 * v * fd.dx + h * D * fd.inner(fd.grad(cp1), fd.grad(v)) * fd.dx
        # hier ist schon die Randbedingung dabei, das j am Rand (3) einen konstanten Strom von 1 hat
        L = c0 * v * fd.dx + fd.Constant(-1.) * h * v * fd.ds(3) # Neumann rbd an surface 3

        bcs = [
                fd.DirichletBC(V, fd.Constant(0.), 4)
        ]


        fd.solve(a == L, c0, bcs) # speichert die Lösung cp1 in c0
        outfile.write(c0) # damit später Animationen gemacht werden können
        #c0 = uu




print(f"successfully written to: {filename}")

#VV = VectorFunctionSpace(mesh, "CG", 1)
#j = project(- sigma * grad(uu), VV)
#VTKFile("j.pvd").write(j)

#VVD = VectorFunctionSpace(mesh, "DG", 0)
#jd = project(- sigma * grad(uu), VVD)
#VTKFile("jd.pvd").write(jd)

