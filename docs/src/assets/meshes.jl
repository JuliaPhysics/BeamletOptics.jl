using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3

include(joinpath(@__DIR__, "render_utils.jl"))

##
benchy = IntersectableObject(joinpath(@__DIR__, "Benchy.stl"))
mesh = BMO.shape(benchy)

translate3d!(benchy, [0, 45mm, -18mm])
zrotate3d!(benchy, deg2rad(41))

s = System([benchy])
pos = [0,0,0]
dir = [0,1,0]
b = Beam(pos, dir)
ray = first(BMO.rays(b))

solve_system!(s, b)

cview = [
  0.443883   0.896085   1.38778e-17  -0.0325398
 -0.183446   0.0908714  0.978821     -0.00727407
  0.877106  -0.434482   0.204719     -0.0461568
  0.0        0.0        0.0           1.0
]

cview = [
  0.244077   0.969756  -1.66533e-16  -0.0394889
 -0.425623   0.107125   0.898537     -0.0106579
  0.871362  -0.219312   0.438897     -0.0577219
  0.0        0.0        0.0           1.0
]

fig = Figure(size=(600,350))
ax = LScene(fig[1,1])

poly!(ax, BMO.vertices(mesh), BMO.faces(mesh), strokewidth=.1, shading=true, transparency=true, alpha=.05, color=:grey)
poly!(ax, BMO.vertices(mesh), BMO.faces(mesh)[251, :]; color=:red)
poly!(ax, BMO.vertices(mesh), BMO.faces(mesh)[37, :]; color=:red, transparency=true, alpha=0.25)
poly!(ax, BMO.vertices(mesh), BMO.faces(mesh)[45, :]; color=:red, transparency=true, alpha=0.25)
poly!(ax, BMO.vertices(mesh), BMO.faces(mesh)[59, :]; color=:red, transparency=true, alpha=0.25)
poly!(ax, BMO.vertices(mesh), BMO.faces(mesh)[123, :]; color=:red, transparency=true, alpha=0.25)
poly!(ax, BMO.vertices(mesh), BMO.faces(mesh)[192, :]; color=:red, transparency=true, alpha=0.25)


render!(ax, b; show_pos=true)
arrow!(ax, pos, dir)

render!(ax, ray; linestyle=:dashdot, flen=200mm)

display(fig)
hide_axis(ax)
set_view(ax, cview)
save("mtalgorithm.png", fig; px_per_unit=4, update = false)