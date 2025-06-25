using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

## collimated source
cs = CollimatedSource([0,0,0], [0,1,0], 15e-3; num_rings=3, num_rays=200)

cs_fig = Figure(; size=(600,200))
display(cs_fig)
cs_ax = LScene(cs_fig[1,1])
hide_axis(cs_ax)

render!(cs_ax, cs, show_pos=true, flen=0.05, color=:red, render_every=1)

cs_view = [
    0.26513     0.964213    7.68482e-16  -0.0180273
    0.0488896  -0.0134432   0.998714      0.00108822
    0.962972   -0.264789   -0.0507042    -0.025119
    0.0         0.0         0.0           1.0
]

set_view(cs_ax, cs_view)
save("collimated_beam_source.png", cs_fig; px_per_unit=8, update = false)

## uniform disc source
cs = UniformDiscSource([0,0,0], [0,1,0], 15e-3; N=100)

cs_fig = Figure(; size=(600,200))
display(cs_fig)
cs_ax = LScene(cs_fig[1,1])
hide_axis(cs_ax)

render!(cs_ax, cs, show_pos=true, flen=0.05, color=:red, render_every=1)

cs_view = [
    0.26513     0.964213    7.68482e-16  -0.0180273
    0.0488896  -0.0134432   0.998714      0.00108822
    0.962972   -0.264789   -0.0507042    -0.025119
    0.0         0.0         0.0           1.0
]

set_view(cs_ax, cs_view)
save("collimated_uniform_beam_source.png", cs_fig; px_per_unit=8, update = false)
