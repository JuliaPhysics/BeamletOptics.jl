using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

##
c_view = [
  0.605459    0.783477  -0.139944  -0.0298959
 -0.0484677   0.211806   0.976109  -0.0106539
  0.7944     -0.584211   0.166213  -1.69474
  0.0         0.0        0.0        1.0
]

gb = GaussianBeamlet([0,0,0], [0,1,0], 9e-4, 2e-3; support=[0,0,1])

l = ThinLens(20e-3, 20e-3, BMO.inch, 1.5)
s = System([l])

translate3d!(l, [0, 30e-3, 0])

fig = Figure(size=(600,300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, gb; flen=0.1-0.03, show_beams=true, show_pos=true, color=RGBAf(1,0,0,0.1))

set_orthographic(ax)
set_view(ax, c_view)


save("gbtest1.png", fig; px_per_unit=8, update = false)

##
solve_system!(s, gb)

fig = Figure(size=(600,300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, gb; flen=0.1-0.03, show_beams=true, show_pos=true, color=RGBAf(1,0,0,0.1))
render!(ax, l; color=RGBAf(1, 1, 1, .1))

c_view = [
  0.653709   0.756746  -5.27356e-16  -0.0173098
 -0.352701   0.304679   0.884744     -0.00702154
  0.669526  -0.578366   0.466077     -0.0223342
  0.0        0.0        0.0           1.0
]

set_view(ax, c_view)
save("gbtest2.png", fig; px_per_unit=8, update = false)

##
using CairoMakie

CairoMakie.activate!()

zs = 0:1e-5:0.125
w, R, ψ, w0 = BeamletOptics.gauss_parameters(gb, zs)

params = Figure(size=(600, 250))
waist = Axis(params[1,1], ylabel="w(z) [m]")
lines!(waist, zs, w, color=:blue)

curv = Axis(params[2,1], yscale=log10, yaxisposition = :right, xlabel="z (along beam optical axis) [m]", ylabel="|R(z)| [m⁻¹]")
lines!(curv, zs, abs.(R) .+ 1, color=:red)

linkxaxes!(waist, curv)
hidexdecorations!(waist, grid=false)
rowgap!(params.layout, 1, 0)

save("gauss_parameters.png", params; px_per_unit=8, update = false)