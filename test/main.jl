using SCDI
using LinearAlgebra
using FileIO, GLMakie

cube = SCDI.Geometry(load("test\\cube.stl"))
SCDI.translate3d!(cube, [0.01,-0.01,-0.01])
SCDI.zrotate3d!(cube, π/4)
SCDI.xrotate3d!(cube, π/4)
SCDI.translate3d!(cube, [0,0.005,0])

f = Figure()
ax = f[1,1] = Axis3(f,aspect = (1,1,1))
xlims!(ax, (-0.025, 0.025))
ylims!(ax, (-0.025, 0.025))
zlims!(ax, (-0.025, 0.025))
mesh!(cube.vertices, cube.faces)
f

ray = SCDI.Ray([-1,0,0], [1,0,0])
t = intersect3d(cube, ray)

lines!([-1,-1+t],[0,0],[0,0], color=:red)
f

mesh!(cube.vertices[3*(2-1)+1:2*3,:]', [1 2 3])
f

# true intersect faces 10 42
# fake intersect faces 2 8

faces = copy(cube.vertices[3*(2-1)+1:2*3,:])
intersect3d(faces, ray, kϵ=1e-6)