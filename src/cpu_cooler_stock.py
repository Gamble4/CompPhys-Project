import firedrake as fd
import math

###
###

def processFunction(t, tmax=0., sigma=1, C=10., offset=10.):
   o = C * math.exp(-1/2 * (t-tmax)**2 / sigma**2) + offset
   return fd.Constant(o)

filename = "cpu_cooler_stock.pvd"
#mesh = fd.Mesh("../blender_salome/cpu_cooler_stock/mesh.msh")
mesh = fd.Mesh("../meshes/bone_with_impurity.msh")

print("setting up function space and functions, ... ", end="")

V = fd.FunctionSpace(mesh, "CG", 1) # testfunktionen, CG = continuous galerkin, Hutbasis

v = fd.TestFunction(V) # 
c = fd.TrialFunction(V) # nach zeitschritt

# use this for c0, set up initial c distribution c(t=0, x)
VD = fd.FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten
c0 = fd.Function(V)

x, y, z = fd.SpatialCoordinate(mesh)
Cw = fd.Constant(20.) # Wassertemp
CHIPT0 = fd.Constant(40.)
alpha = fd.Constant(4) # Koeffizient der Konvektion, kann auch Ortsabh. sein

# setze Anfangstemperatur
c0.interpolate(Cw, subset=mesh.cell_subset(1)) 
c0.interpolate(Cw, subset=mesh.cell_subset(2)) 

# set Diffusionskoeffizient 
D = fd.Function(VD)
D.interpolate(fd.Constant(0.5),  subset=mesh.cell_subset(1))
D.interpolate(fd.Constant(1.0),  subset=mesh.cell_subset(2))

outfile = fd.VTKFile(filename)
print("done")


# chip BCS
bcs = [                                
        fd.DirichletBC(V, fd.Constant(CHIPT0),3),
        ]
bc = fd.DirichletBC(V, fd.Constant(CHIPT0),3)
# wir setzen adiabatische Randbedingungen, bis auf einen konvektiven Rand
t = 0.
dt = 0.1
h = fd.Constant(dt) # Zeitschritt für implizite Euler Integration
dx = fd.Measure("dx", domain=mesh)


q = fd.Function(V).assign(c0)
q_g = fd.Function(V)
bc.apply(q)
q_g.interpolate(Cw)
bc.apply(q_g)
q2 = fd.Function(V)
dq = fd.Function(V)

a = fd.inner(c, v) * dx
L = -h * D * fd.inner(fd.grad(q), fd.grad(v)) * dx #
    #- h * D * fd.inner(fd.grad(q_g), fd.grad(v)) * dx \


prob1 = fd.LinearVariationalProblem(a, L, dq)
solv1 = fd.LinearVariationalSolver(
                prob1,
                solver_parameters={
                    "ksp_type" : "cg",
                    "pc_type" : "jacobi"
                    })

for i in range(10):

        
    print(f"iteration: {i}", end="\r")
    solv1.solve()
    q.assign(q + dq)



    outfile.write(q) # damit später Animationen gemacht werden können
    t += dt


print(f"successfully written to: {filename}")

