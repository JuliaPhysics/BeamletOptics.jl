using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

const mm = 1e-3

##
fig = Figure(size=(600,300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

# coord. sys
arrow!(ax, [0,0,0], [1,0,0]; color=:red)
arrow!(ax, [0,0,0], [0,1,0]; color=:green)
arrow!(ax, [0,0,0], [0,0,1]; color=:yellow)
text!(Point3(5mm, 0, 1mm), text="x")
text!(Point3(0, 5mm, 1mm), text="y")
text!(Point3(0, 0mm, 6.5mm), text="z")

# ray
cube = IntersectableObject(BMO.CubeMesh(5mm))
translate_to3d!(cube, [-20mm, 20mm, 10mm])
xrotate3d!(cube, deg2rad(30))
yrotate3d!(cube, deg2rad(30))
zrotate3d!(cube, deg2rad(-30))

ray_pos = [0mm,5mm,5mm]
ray_pos = [5mm,1mm,10mm]
ray_dir = position(cube) - ray_pos + [5mm, 0, 1mm]
ray = Ray(ray_pos, ray_dir)

ray_dir = BMO.direction(ray)

is = BMO.intersect3d(cube, ray)
hit_pos = ray_pos + length(is) * ray_dir

arrow!(ax, ray_pos, ray_dir; scale=0.5)
scatter!(ax, Point3(hit_pos); color=:black)

scatter!(ax, Point3(ray_pos); color=:blue)

arrow!(ax, hit_pos, BMO.normal3d(is); scale=0.5, color=:black)

render!(ax, ray; flen=length(is), linestyle=:dashdot)
render!(ax, cube; color=lens_color(0.5))

c_view = [
 -0.685229    0.728328   6.93889e-17  -0.012345
 -0.0391483  -0.0368317  0.998554     -0.00716188
  0.727275    0.684239   0.053751     -0.0278491
  0.0         0.0        0.0           1.0
]

set_view(ax, c_view)
save("ray_showcase.png", fig; px_per_unit=8, update = false)