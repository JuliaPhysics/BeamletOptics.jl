using GLMakie, BeamletOptics

const mm = 1e-3

##
lens = ThinLens(50mm, 50mm, BeamletOptics.inch, 1.5)

sd = Spotdetector(10e-3)
translate3d!(sd, [0, 50mm, 0])

system = System([lens, sd])

## render system
system_fig = Figure(size=(600,300))
limits = (-0.02, 0.02, -0.025, 0.05, -0.015, 0.015)
aspect = (1,1.875,0.75)
system_ax = Axis3(system_fig[1,1], aspect=aspect, limits=limits, azimuth=-1, elevation=0.)

hidedecorations!(system_ax)
hidespines!(system_ax)

render!(system_ax, system)

beam = Beam([0,-50mm,0], [0,1,0], 1e-6)

aperture = BeamletOptics.inch*0.8

zs = LinRange(-aperture/2, aperture/2, 25)

for z in zs
    x = BeamletOptics.position(first(beam.rays))[1]
    y = BeamletOptics.position(first(beam.rays))[2]
    BeamletOptics.position!(first(beam.rays), Point3{Float64}(x, y, z))
    solve_system!(system, beam)
    render!(system_ax, beam, show_pos=true)
end

## render diagram
n_rings = 20
n_rays = 5000
BeamletOptics.create_spot_diagram(system, beam, aperture; n_rings, n_rays)
# Rerun for time without compile
t1 = @timed BeamletOptics.create_spot_diagram(system, beam, aperture; n_rings, n_rays)

spot_fig = Figure(size=(600,400))
spot_ax = Axis(spot_fig[1,1], aspect=1, xlabel="x [mm]", ylabel="y [mm]")
sc = scatter!(spot_ax, sd.data, markersize=3, color=:blue)

extime = trunc(t1.time*1e3, digits=2)

leg_string = "
   # of traces: $n_rays \n
   # of rings: $n_rings \n
   Ex. Time: $extime ms \n
"

Legend(spot_fig[1,2], [sc], [leg_string], "Quick stats.")