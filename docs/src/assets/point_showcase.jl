using GLMakie, BeamletOptics

const BMO = BeamletOptics

GLMakie.activate!(; ssao=true)

Base.include(@__MODULE__, "render_utils.jl")

## point source
ps = PointSource([0,0,0], [0,1,0], deg2rad(10); num_rings=3, num_rays=200)

ps_fig = Figure(; size=(600,200))
display(ps_fig)
ps_ax = LScene(ps_fig[1,1])

for i = 1:1:length(BMO.beams(ps))
    render_beam!(ps_ax, BMO.beams(ps)[i], show_pos=true, flen=0.1, color=:red)
end

ps_view = [
    0.00272773   0.999996     5.03287e-16  -0.0584386
    0.424361    -0.00115755   0.905492      2.95341e-5
    0.905489    -0.00246994  -0.424362     -0.04708
    0.0          0.0          0.0           1.0
]

set_view(ps_ax, ps_view)

hide_axis(ps_ax)