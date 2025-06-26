using GLMakie, BeamletOptics

pd = Photodetector(1e-3, 1000)
pd_body = MeshDummy(joinpath(asset_dir, "FDS010.stl"))
zrotate3d!(pd_body, π)
translate3d!(pd, [0,-3e-3,0])

pd_group = ObjectGroup([pd, pd_body])

system = System([pd])

g1 = GaussianBeamlet([0,-20e-3,0], [0,1,1e-3], 532e-9, 1e-5)
g2 = GaussianBeamlet([0,-20e-3,0], [0,1,0], 532e-9, 5e-4)

ϕ = 0
g2.E0 *= exp(im*ϕ)

solve_system!(system, g1)
solve_system!(system, g2)

## render fringes
fringes_fig = Figure()
heat = Axis(fringes_fig[1, 1], xlabel="x [mm]", ylabel="y [mm]", aspect=1)
hm = heatmap!(heat, pd.x*1e3, pd.y*1e3, intensity(pd), colormap=:viridis)
cb = Colorbar(fringes_fig[1, 2], hm, label="Intensity [W/m²]")

## render system
detector_fig = Figure(size=(600, 280))
aspect = (.5,1.5,.5)
limits = (-0.005, 0.005, -0.02, 0.01, -0.005, 0.005)
rend = Axis3(detector_fig[1,1], aspect=:data, limits=limits, azimuth=5.43, elevation=0, perspectiveness=0.1)

hidedecorations!(rend)
hidespines!(rend)

render!(rend, pd_body)
render!(rend, pd)
render!(rend, g2)
render!(rend, g1)
