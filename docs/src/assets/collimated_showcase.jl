using GLMakie, BeamletOptics

const BMO = BeamletOptics

GLMakie.activate!(; ssao=true)

Base.include(@__MODULE__, "render_utils.jl")

## collimated source
cs = CollimatedSource([0,0,0], [0,1,0], 15e-3; num_rings=3, num_rays=200)

cs_fig = Figure(; size=(600,200))
display(cs_fig)
cs_ax = LScene(cs_fig[1,1])

for i = 1:1:length(BMO.beams(cs))
    render_beam!(cs_ax, BMO.beams(cs)[i], show_pos=true, flen=0.05, color=:red)
end

cs_view = [
    0.26513     0.964213    7.68482e-16  -0.0180273
    0.0488896  -0.0134432   0.998714      0.00108822
    0.962972   -0.264789   -0.0507042    -0.025119
    0.0         0.0         0.0           1.0
]

set_view(cs_ax, cs_view)

hide_axis(cs_ax)