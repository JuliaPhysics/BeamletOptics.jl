using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

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
num_rings = 20
num_rays = 5000

cs = CollimatedSource([0,-50mm,0], [0,1,0], aperture, 1e-6; num_rings, num_rays)

t1 = @timed solve_system!(system, cs)

render!(system_ax, cs, color=:blue, show_pos=false, render_every=50)

save("spot_diagram_system.png", system_fig, px_per_unit=4)

## render diagram
spot_fig = Figure(size=(600,400))
spot_ax = Axis(spot_fig[1,1], aspect=1, xlabel="x [mm]", ylabel="y [mm]")
sc = scatter!(spot_ax, sd.data, markersize=3, color=:blue)

extime = trunc(t1.time*1e3, digits=2)

leg_string = "
   # of traces: $num_rays \n
   # of rings: $num_rings \n
   Ex. Time: $extime ms \n
"

Legend(spot_fig[1,2], [sc], [leg_string], "Quick stats.")

save("spot_diagram_showcase.png", spot_fig, px_per_unit=4)