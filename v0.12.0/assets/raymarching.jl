using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3

include(joinpath(@__DIR__, "render_utils.jl"))

##
box_dim = 8mm
box = BMO.BoxSDF(box_dim, box_dim, box_dim)
ring = BMO.RingSDF(15mm, 2mm, 2mm)
cylinder = BMO.CylinderSDF(5mm, 5mm)
sphere = BMO.SphereSDF(10mm)

translate3d!(box,       [-10mm, 20mm, 10mm])
translate3d!(ring,      [0, 50mm, 0])
translate3d!(cylinder,  [15mm, 35mm, -5mm])
translate3d!(sphere,    [-12mm, 80mm, 0])

xrotate3d!(cylinder, deg2rad(90))
yrotate3d!(cylinder, deg2rad(-20))

xrotate3d!(box, deg2rad(30))
yrotate3d!(box, deg2rad(30))

cview = [
 -0.431147   0.902282  8.32667e-17  -0.0544151
 -0.420532  -0.200948  0.884744      0.0126278
  0.798289   0.381455  0.466077     -0.108933
  0.0        0.0       0.0           1.0
]

##
fig = Figure(size=(600,280))
ax = LScene(fig[1,1])

render!(ax, box)
render!(ax, ring)
render!(ax, cylinder)
render!(ax, sphere)

##
pos = [0,-10mm,0]
dir = [0,1,0]

sdfs = [box, ring, cylinder, sphere]

global new_pos = pos

max_i = 16
for i = 1:max_i
    distances = BMO.sdf.(sdfs, Ref(new_pos))
    step = minimum(distances)
    if i != max_i
        render_sphere!(ax, new_pos, step, color=RGBAf(1,1,1,0.25))
        ray = Ray(new_pos, dir)
        render!(ax, ray, show_pos=true, flen=step)
    else
        arrow!(ax, new_pos, dir; color=:blue)
        scatter!(ax, Point3(new_pos); color=:blue)
    end
    global new_pos = new_pos .+ dir * step
end

display(fig)
hide_axis(ax)
set_view(ax, cview)
save("raymarching.png", fig; px_per_unit=4, update = false)
