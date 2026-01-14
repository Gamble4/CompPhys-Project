from firedrake import *
from petsc4py import PETSc


# ----------------------------
# 4. PETSc RHS-Funktion
#    M u_t = -K u
# ----------------------------
def rhs(ts, t, x, f):
    """
    f = F(u) = -K u
    mult performs product K@x -> f gespeichert
    """
    f.zeroEntries()
    K.mult(x, rhs_vec)
    rhs_vec.scale(-1.0)
    rhs_vec.axpy(1.0, b)

    # löse M f = rhs_vec
    kspM.solve(rhs_vec, f)
# ----------------------------
# 5. Jacobian
#    J = a*M - K
# ----------------------------
def jac(ts, t,  a, J, P):
    """
    J = dF/du + a*M
    """
    J.zeroEntries()
    J.axpy(1, M, structure=PETSc.Mat.Structure.SAME)
    J.axpy(-1.0, K, structure=PETSc.Mat.Structure.SAME)
    J.assemble()
# ----------------------------
# 1. Mesh und Funktionsraum
# ----------------------------
mesh = Mesh("mesh.msh")
V = FunctionSpace(mesh, "CG", 1)

# Lösung
u0 = Function(V)
u = TrialFunction(V)
v = TestFunction(V)
x, y, z = SpatialCoordinate(mesh)
u0.interpolate(sin(pi * x) * sin(pi * y), subset=mesh.cell_subset(1))

bc1 = DirichletBC(V, Constant(10.), 2)
alpha = Constant(5.)
Cw = Constant(0.)
VD = FunctionSpace(mesh, "DG", 0) # DG = discrete Galerkin, für konstanten

D = Function(VD)
D.interpolate(Constant(1.), subset=mesh.cell_subset(1))

# Massenmatrix
M_form = inner(u, v) * dx

# Steifigkeitsmatrix (Diffusion)
K_form = D * inner(grad(u), grad(v)) * dx + alpha * u * v *ds(3)
b_form = alpha * Cw * v * ds(3)


# Assemblieren als PETSc-Matrizen
M = assemble(M_form, mat_type="aij").petscmat
K = assemble(K_form, mat_type="aij").petscmat

kspM = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
kspM.setOperators(M)
kspM.setType("cg")
kspM.getPC().setType("jacobi")   # M ist SPD 
kspM.setUp()


with assemble(b_form).dat.vec as b:

    rhs_vec = M.createVecRight()
    def rhs(ts, t, x, f):
        """
        f = F(u) = -K u
        mult performs product K@x -> f gespeichert
        """
        f.zeroEntries()
        K.mult(x, rhs_vec)
        rhs_vec.scale(-1.0)
        rhs_vec.axpy(1.0, b)

        # löse M f = rhs_vec
        kspM.solve(rhs_vec, f)
    # ---------------------------


    ts = PETSc.TS().create(comm=PETSc.COMM_WORLD)

    ts.setProblemType(PETSc.TS.ProblemType.LINEAR)
    ts.setType(PETSc.TS.Type.RK)

    outfile = VTKFile("diff_lin.pvd")

    # Zeitparameter
    ts.setTime(0.0)
    ts.setTimeStep(1e-1)
    ts.setMaxTime(10)
    ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
    apply_bc = lambda ts, bc : bc1.apply(u0)
    ts.setPostStep(apply_bc, bc1)
    print_t = lambda ts : print(ts.getTime(), end="\r")
    ts.setPostStep(print_t)
    
    ts.setRHSFunction(rhs)
    

    print("solving")
    with u0.dat.vec as x: 
        ts.solve(x)

    outfile.write(u0)
    PETSc.Sys.Print("Simulation beendet.")


