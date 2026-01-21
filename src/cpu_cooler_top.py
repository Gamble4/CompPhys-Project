import firedrake as fd
import time


def pF(t, tmax, sigma, C, offset=0.):
    return C * fd.exp(-1/2 * (t-tmax)**2 / sigma**2) + offset
    

filename = "../solutions/test.pvd"
mesh = fd.Mesh("../../CompPhys-Proj/blender_salome/cpu_cooler_simple/mesh.msh")
x, y, z = fd.SpatialCoordinate(mesh)

V = fd.FunctionSpace(mesh, "CG", 1)
VD = fd.FunctionSpace(mesh, "DG", 0)

v = fd.TestFunction(V)
u = fd.TrialFunction(V)
c = fd.Function(V) # Temperaturfeld


DPIPE = 5.
DFINS = 0.01
DBLOC = 0.02

C0 = fd.Constant(35.)
CAIR = fd.Constant(20.)
ALPHA = fd.Constant(0.01)

t = 0.
dt = 0.05
T = 20.
step = 0

fd_t = fd.Constant(0.)


D = fd.Function(VD)
D.interpolate(DFINS, subset=mesh.cell_subset(1)) # Fins
D.interpolate(DPIPE, subset=mesh.cell_subset(2)) # pipes
D.interpolate(DBLOC, subset=mesh.cell_subset(3))

c.interpolate(C0)

alpha = fd.Function(VD)
alpha.interpolate(ALPHA, subset=mesh.cell_subset(1)) # Fins
alpha.interpolate(0., subset=mesh.cell_subset(2))
alpha.interpolate(0., subset=mesh.cell_subset(3))


# Überlagerung an Gaußprofilen für die prozessfunktionen
def source(t):
    return pF(t, 15, 4, 3, 0) + pF(t, 5, 3, 4, 0)

a = (
    u*v*fd.dx
    +  D * dt/2 * fd.inner(fd.grad(u), fd.grad(v)) * fd.dx
    +  dt/2 * alpha * u * v * fd.ds(5)
    )

L = (
        c*v*fd.dx
        - D * dt/2 * fd.inner(fd.grad(c), fd.grad(v)) * fd.dx
        - dt/2 * alpha * c * v * fd.ds(5)
        + dt * alpha * CAIR * v * fd.ds(5)
        + dt/2 * ( source(fd_t) + source(fd_t + dt) ) * v * fd.dx(3) # Quellterm
    )

sparams = {
    "ksp_type": "cg",
    "pc_type": "hypre",
    "ksp_monitor" : None,
    "ksp_view" : None,
    "pc_view" : None
}

# für bcs : bcs=bcs in LinearVariationalSolver einfügen
problem = fd.LinearVariationalProblem(a, L, c)
solver = fd.LinearVariationalSolver(
    problem,
    solver_parameters=sparams,
    )

outfile = fd.VTKFile(filename)

start = time.time()
while t <= T:

    solver.solve()
    fd_t.assign(t)

    print(f"time: {t:e} : {t/T * 100:.2f}%", end="\r")
    t += dt
    step += 1
    if step % 10 == 0:
        outfile.write(c)

#outfile.write(c)

print(f"took {time.time() - start:.2f}s to finish")
print(f"done writing to {filename}")

