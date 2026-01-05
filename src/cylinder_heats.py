from firedrake import *

mesh = Mesh("../meshes/cylinder_heats.msh")

V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

VD = FunctionSpace(mesh, "DG", 0)
sigma = Function(VD)
sigma.interpolate(Constant(1.), subset=mesh.cell_subset(1))


a = sigma * inner(grad(u), grad(v)) * dx
L = Constant(1.) * v * ds(3)

bcs = [
        DirichletBC(V, Constant(10.), 2)
]

uu = Function(V)

solve(a == L, uu, bcs)
VTKFile("u.pvd").write(uu)

VV = VectorFunctionSpace(mesh, "CG", 1)
j = project(- sigma * grad(uu), VV)
VTKFile("j.pvd").write(j)

VVD = VectorFunctionSpace(mesh, "DG", 0)
jd = project(- sigma * grad(uu), VVD)
VTKFile("jd.pvd").write(jd)
