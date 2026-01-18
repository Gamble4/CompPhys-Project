import firedrake as fd
import math

###
# tags

###

def pF(t, tmax=0., sigma=1, C=10., offset=0.):
   return C * math.exp(-1/2 * (t-tmax)**2 / sigma**2) + offset
   

filename = "ProtoChipCooled.pvd"
mesh = fd.Mesh("../meshes/ProtoChipCooled.msh")


print("setting up function space and functions, ... ", end="")

V = fd.FunctionSpace(mesh, "CG", 1) # testfunktionen, CG = continuous galerkin, Hutbasis

v = fd.TestFunction(V) # 
u = fd.TrialFunction(V) # nach zeitschritt

# use this for c0, set up initial c distribution c(t=0, x)
VD = fd.FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten
c = fd.Function(V)

x, y, z = fd.SpatialCoordinate(mesh)
Cw = fd.Constant(20.) # Wassertemp
CHIPT0 = fd.Constant(40.)
alpha = fd.Constant(0.1) # Koeffizient der Konvektion, kann auch Ortsabh. sein

# setze Anfangstemperatur
c.interpolate(CHIPT0)


# set Diffusionskoeffizient 
D = fd.Function(VD)
D.interpolate(fd.Constant(0.01),  subset=mesh.cell_subset(1))
D.interpolate(fd.Constant(0.02),  subset=mesh.cell_subset(2))
D.interpolate(fd.Constant(0.02),  subset=mesh.cell_subset(3))
D.interpolate(fd.Constant(0.02),  subset=mesh.cell_subset(4))
D.interpolate(fd.Constant(1.),  subset=mesh.cell_subset(5))
D.interpolate(fd.Constant(1.),  subset=mesh.cell_subset(6))
D.interpolate(fd.Constant(0.01),  subset=mesh.cell_subset(7))

outfile = fd.VTKFile(filename)
print("done")

fd_t = fd.Constant(0.)

bcval = pF(fd_t, 10, 1, 10, CHIPT0) + pF(fd_t, 5, 1, 20, CHIPT0)

# chip BCS
bcs = [                                
        fd.DirichletBC(V, bcval, [9, 10, 11])
        ]
# wir setzen adiabatische Randbedingungen, bis auf einen konvektiven Rand
T = 15.
t = 0.
dt = 0.01
h = fd.Constant(dt) # Zeitschritt für implizite Euler Integration
dx = fd.Measure("dx", domain=mesh)

step = 0



a = (
    u*v*fd.dx
    +  dt/2 * fd.inner(fd.grad(u), fd.grad(v)) * fd.dx
    +  dt/2 * alpha * u * v * fd.ds(8)
    )

L = (
        c*v*fd.dx
        - dt/2 * fd.inner(fd.grad(c), fd.grad(v)) * fd.dx
        - dt/2 * alpha * c * v * fd.ds(8)
        + dt * alpha * CHIPT0 * v * fd.ds(8)
        #+ dt/2 * ( pF(fd_t, 10, 1, 10) + pF(fd_t+dt, 10, 1, 10) ) * v * fd.dx # Quellterm
    )

sparams = {
    "ksp_type": "cg",
    "pc_type": "hypre"
}

# für bcs : bcs=bcs in LinearVariationalSolver einfügen
problem = fd.LinearVariationalProblem(a, L, c, bcs=bcs)
solver = fd.LinearVariationalSolver(
    problem,
    solver_parameters=sparams,
    )


while t < T:

    solver.solve()
    fd_t.assign(t)

    print(f"time: {t:e} : {t/T * 100:.2f}%", end="\r")
    t += dt
    step += 1
    if step % 10 == 0:
        outfile.write(c)



print(f"successfully written to: {filename}")

#VV = VectorFunctionSpace(mesh, "CG", 1)
#j = project(- sigma * grad(uu), VV)
#VTKFile("j.pvd").write(j)

#VVD = VectorFunctionSpace(mesh, "DG", 0)
#jd = project(- sigma * grad(uu), VVD)
#VTKFile("jd.pvd").write(jd)

