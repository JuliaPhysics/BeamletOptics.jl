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
filter = PolarizationFilter(20mm)
translate_to3d!(filter, position(m2) + [0, 70mm, 0])
yrotate3d!(filter, deg2rad(0))

##
system = System([m1, m2, m3, filter])

E_circ = [0, 2im, -1] / sqrt(2)
E_lin = [0,1,0]

ray = PolarizedRay(position(m1) + [80mm, 0, 0], [-1, 0, 0], 1000nm, E_lin)
beam = Beam(ray)

solve_system!(system, beam)

##
fig = Figure()#size=(600,300))
ax = LScene(fig[1,1], show_axis=false)

render!(ax, periscope; transparency=true, alpha=0.05)

render!(ax, m1; color=:gold)
render!(ax, m2; color=:gold)
render!(ax, m3; color=:gold)

render!(ax, filter)
# render_lcs!(ax, filter; scale=3)

render!(ax, beam; flen=80mm, color=:blue, show_polarization=true, pol_λ=5mm, pol_color=:blue, pol_amplitude=10mm)

fig