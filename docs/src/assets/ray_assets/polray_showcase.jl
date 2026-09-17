using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

PSF10_03_P01() = RoundPlanoMirror(25.4mm, 6mm)

##
periscope = MeshDummy(joinpath(@__DIR__, "periscope.stl"))

##
m1 = PSF10_03_P01()
m2 = PSF10_03_P01()
m3 = PSF10_03_P01()

translate_to3d!(m1, [24.332mm, 24.13mm, 24.332mm])
zrotate3d!(m1, deg2rad(90))
yrotate3d!(m1, deg2rad(-45))

translate_to3d!(m2, position(m1) + [0, 0, 127.856mm])
xrotate3d!(m2, deg2rad(45+90))

translate_to3d!(m3, position(m2) + [0, 127.856mm, 0])
zrotate3d!(m3, deg2rad(-45))

##
flt = PolarizationFilter(20mm)
translate_to3d!(flt, position(m2) + [0, 70mm, 0])
yrotate3d!(flt, deg2rad(0))

##
system = System([m1, m2, m3, flt])

E_circ = [0, 1im, 1] / sqrt(2)
E_lin = [0,1,0]

start_pos = position(m1) + [150mm, 0, 0]

ray = PolarizedRay(start_pos, [-1, 0, 0], 1000nm, E_lin)
beam = Beam(ray)

beam = AstigmaticGaussianBeamlet(start_pos, [-1, 0, 0], 1000nm, 4mm, E0=E_circ)

solve_system!(system, beam)

##
cview = [
 -0.59507    0.803674  9.15934e-16  -0.0651614
 -0.548087  -0.405824  0.731373     -0.02743
  0.587786   0.435218  0.681977     -0.384558
  0.0        0.0       0.0           1.0
]

fig = Figure(size=(600,450))
display(fig)
ax = LScene(fig[1,1], show_axis=false)

render!(ax, periscope; transparency=true, alpha=0.05)

render!(ax, m1; color=:gold)
render!(ax, m2; color=:gold)
render!(ax, m3; color=:gold)

render!(ax, flt)
# render_lcs!(ax, flt; scale=3)

render!(ax, beam; flen=200mm, color=RGBAf(1,0,0,0.25), show_polarization=true, pol_λ=5mm, pol_color=:red)

set_view(ax, cview)

save("polray_showcase.png", fig; px_per_unit=4, update = false)