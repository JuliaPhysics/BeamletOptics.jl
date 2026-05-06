using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

const mm = 1e-3

include(joinpath(@__DIR__, "..", "render_utils.jl"))

##
c_view = [
  0.605459    0.783477  -0.139944  -0.0298959
 -0.0484677   0.211806   0.976109  -0.0106539
  0.7944     -0.584211   0.166213  -1.69474
  0.0         0.0        0.0        1.0
]

agb = AstigmaticGaussianBeamlet([0,0,0], [0,1,0], 9e-4, 2e-3; support=[0,0,1])

l = ThinLens(20e-3, 20e-3, BMO.inch, 1.5)
s = System([l])

translate3d!(l, [0, 30e-3, 0])

fig1 = Figure(size=(600,300))
display(fig1)
ax = LScene(fig1[1,1])
hide_axis(ax)

render!(ax, agb; flen=0.1-0.03, show_beams=true, show_pos=true, color=RGBAf(1,0,0,0.1))

set_orthographic(ax)
set_view(ax, c_view)

save("agbtest1.png", fig1; px_per_unit=8, update = false)

##
solve_system!(s, agb; check_invariant=false)

fig = Figure(size=(600,300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, agb; flen=0.1-0.03, show_beams=true, show_pos=true, color=RGBAf(1,0,0,0.1))
render!(ax, l; color=RGBAf(1, 1, 1, .1))

c_view = [
  0.653709   0.756746  -5.27356e-16  -0.0173098
 -0.352701   0.304679   0.884744     -0.00702154
  0.669526  -0.578366   0.466077     -0.0223342
  0.0        0.0        0.0           1.0
]

set_view(ax, c_view)
save("agbtest2.png", fig; px_per_unit=8, update = false)

##
s1 = BMO.CylinderSDF(BMO.inch/2, BMO.inch)
s2 = BMO.CylinderSDF(BMO.inch/2, BMO.inch)
l1 = BMO.Lens(s1, n->1.5)
l2 = BMO.Lens(s2, n->1.5)
translate3d!(l1, [0, 50mm,  0])
translate3d!(l2, [0, 100mm, 0])
xrotate3d!(l1, deg2rad(90))
xrotate3d!(l2, deg2rad(90))
yrotate3d!(l1, deg2rad(90))
translate3d!(l2, [5mm, 0, 0])

system = System([l1, l2])

z_os = 0e-3
w0 = 2mm
λ0 = 1000e-9
dir = [0,1,0]
E0 = [0,0,1]
beam = AstigmaticGaussianBeamlet([0., 0, z_os], dir, λ0, w0; E0, support=[0,0,1])

solve_system!(system, beam; check_invariant=false)

##
fig = Figure(size=(600,300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, l1; color=RGBAf(1, 1, 1, .2))
render!(ax, l2; color=RGBAf(1, 1, 1, .2))
render!(ax, beam; color=RGBAf(1,0,0,0.2), z_res=4, r_res=20, show_waist=true, markersize=4, flen=0.1)

view = [
  0.714009   0.700137  -1.38778e-16  -0.0666259
 -0.219833   0.224189   0.949428     -0.0103548
  0.664729  -0.6779     0.313986     -0.0439393
  0.0        0.0        0.0           1.0
]

set_view(ax, view)
save("agbtest3.png", fig; px_per_unit=8, update = false)