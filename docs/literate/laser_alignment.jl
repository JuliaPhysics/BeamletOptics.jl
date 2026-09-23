# # Laser alignment
#
# This tutorial walks through a small lab-style task: lifting a HeNe laser beam onto a table with a periscope, discovering and correcting a mirror mounting error with an alignment card, and finally focusing the beam onto a camera sensor.
#
# ```@raw html
# <div class="bmo-card">
# <p><Badge type="tip" text="Beginner" /></p>
# ```
#
# You will learn how to:
#
# 1. Place components with the kinematic API
# 2. Trace a [`Beam`](@ref) and read out a [`Detector`](@ref)
# 3. Find and correct a misaligned mirror
# 4. Focus a [`GaussianBeamlet`](@ref) and measure its waist
#
# ```@raw html
# </div>
# ```
#
# ```@raw html
# <div class="bmo-card">
# <p><a href="../laser_alignment.jl">Download tutorial <code>.jl</code></a> &middot; <a href="../laser_alignment.ipynb">Download Jupyter notebook</a></p>
# </div>
# ```
#
# Our HeNe laser sits at a beam height of 50 mm, while the rest of the optical table works at 150 mm. A periscope built from two mirrors lifts the beam onto the table height. The first mirror, `M1`, turns out to have a small mounting error; we will find it with an alignment card and correct it. Finally, a lens focuses the beam onto a camera.
#
# !!! info "Units"
#     Unless stated otherwise, this package assumes SI units for input parameters. We define `const mm = 1e-3` below and use `mm` throughout to make lengths easier to read, e.g. `50mm` is 50 millimeters expressed in meters.
#
# !!! tip "Before you start"
#     This tutorial assumes basic familiarity with Julia. If you are new to the language, the [Getting started](https://docs.julialang.org/en/v1/manual/getting-started/) section of the Julia manual and the [Julia learning resources](https://julialang.org/learning/) are good starting points; [Modern Julia Workflows](https://modernjuliaworkflows.org/) covers package environments and editor setup.
#
#     We also assume that you run the code in [VS Code](https://code.visualstudio.com/) with the [Julia extension](https://www.julia-vscode.org/docs/stable/gettingstarted/) installed. Execute the code blocks one after another in the integrated Julia REPL, e.g. by selecting them and pressing `Shift+Enter`. You need to install the packages once beforehand: press `]` in the REPL and type `add BeamletOptics GLMakie`.
#
#     Figures appear only when the figure object is returned or displayed. On this page, the plots are shown as images below the code blocks. When running the code yourself, end each plotting block with the figure variable (e.g. `fig`) or call `display(fig)`. With GLMakie, `display(fig)` opens an interactive window where the 3D scene can be rotated and zoomed. GLMakie handles the 2D plots in this tutorial as well.
#
# ## Building the periscope
#
# We start by defining the two mirrors of the periscope, a [`RoundPlanoMirror`](@ref) with an outer diameter of one inch and a thickness of 6 mm each.

using GLMakie, BeamletOptics
const BMO = BeamletOptics
const mm = 1e-3

h1, h2 = 50mm, 150mm     # beam height of the laser and of the table optics
λ = 632.8e-9             # HeNe wavelength

m1 = RoundPlanoMirror(BMO.inch, 6mm)
m2 = RoundPlanoMirror(BMO.inch, 6mm)
xrotate3d!(m1, deg2rad(-45)); translate3d!(m1, [0, 100mm, h1])
xrotate3d!(m2, deg2rad(135)); translate3d!(m2, [0, 100mm, h2])

# Every component spawns at the global origin facing the +y-axis, which is the direction of the optical axis in this tutorial. Rotating `M1` by −45° around the x-axis tilts its face so that the incoming horizontal beam is sent straight up. `M2` is rotated by 135° rather than −45°, so that its **front** face points back down at the beam coming up from `M1`. With −45° the beam would instead hit the *back* of the 6 mm thick mirror substrate and leave the periscope roughly 10 mm too low, missing the rest of the setup entirely.
#
# We can now trace a single [`Beam`](@ref) through the periscope and render the result. `System`s bundle all components that take part in a simulation, and [`solve_system!`](@ref) performs the actual ray tracing.

system = System([m1, m2])
beam = Beam(Ray([0, 0, h1], [0, 1.0, 0], λ))
solve_system!(system, beam)

fig = Figure(size=(600, 400))
ax = Axis3(fig[1,1], aspect=:data, azimuth=0.3π, elevation=0.15π)
hidedecorations!(ax)
render!(ax, system)
render!(ax, beam, color=:red, flen=0.3)
save("periscope.png", fig, px_per_unit=4); nothing #hide #md
fig #!md

# ![Aligned periscope lifting the beam from 50 mm to 150 mm](periscope.png)
#
# ## Finding the misalignment
#
# In practice, mirror mounts are never perfectly aligned. We simulate a small mounting error on `M1` and place an alignment card 300 mm downstream of it (200 mm past `M2`) to see where the beam actually lands.

xrotate3d!(m1, deg2rad(0.5))           # M1 was mounted 0.5° off

card = Detector(BMO.inch, false)       # alignment card: records hits, lets the beam pass
translate3d!(card, [0, 300mm, h2])

system = System([m1, m2, card])
beam = Beam(Ray([0, 0, h1], [0, 1.0, 0], λ))
empty!(card)
solve_system!(system, beam)
offset = spot_diagram(card)[1]         # [x, z] on the card
println("offset on card: ", round.(offset ./ mm, digits=2), " mm")

# A [`Detector`](@ref) with `stop = false` behaves like a real alignment card: it records where the beam hits, but lets the beam continue propagating through the rest of the system rather than absorbing it. Detectors accumulate hits over successive calls to `solve_system!`, so we call [`empty!`](@ref) beforehand to discard any previous data.

fig2 = Figure(size=(750, 400))
ax2 = Axis3(fig2[1,1], aspect=:data, azimuth=0.3π, elevation=0.15π)
hidedecorations!(ax2)
render!(ax2, system)
render!(ax2, beam, color=:red, flen=0.3)

spot_ax = Axis(fig2[1,2], aspect=1, xlabel="x [mm]", ylabel="z [mm]",
    limits=(-12.7, 12.7, -12.7, 12.7), title="Card")
scatter!(spot_ax, [0.0], [0.0], color=:black, marker=:xcross, markersize=16)
pts = spot_diagram(card)
scatter!(spot_ax, [p[1]/mm for p in pts], [p[2]/mm for p in pts], color=:red)

save("card_misaligned.png", fig2, px_per_unit=4); nothing #hide #md
fig2 #!md

# ![Misaligned beam missing the target on the alignment card](card_misaligned.png)
#
# The black cross marks the intended target at the center of the card, while the red dot shows where the beam actually lands, offset by roughly 5 mm.
#
# ## Correcting the mirror
#
# Tilting a mirror by a small angle ``\delta`` deflects the reflected beam by ``2\delta``. Since the card sits a path length ``L`` behind `M1` (up to `M2`, then across to the card), the resulting offset on the card is approximately
#
# ```math
# \text{offset} = 2\delta L
# ```
#
# We can invert this relation to compute the correction angle from the measured offset and apply it to `M1`.

L = (h2 - h1) + (300mm - 100mm)        # M1 → M2 → card
δ = offset[2] / (2L)
println("correction: ", round(rad2deg(δ), digits=3), "°")
xrotate3d!(m1, δ)

empty!(card)
beam = Beam(Ray([0, 0, h1], [0, 1.0, 0], λ))
solve_system!(system, beam)
println("offset after correction: ", round.(spot_diagram(card)[1] ./ 1e-6, digits=2), " µm")

# The correction angle exactly cancels the 0.5° mounting error we introduced above, and the residual offset on the card drops from several millimeters to well under a micrometer.
#
# ## Focusing onto a camera
#
# With `M1` corrected, we add a plano-convex lens to focus the beam and place a camera at the resulting waist. The lens is modeled with a [`SellmeierEquation`](@ref) dispersion model for N-BK7 glass and built with the [`SphericalLens`](@ref) convenience constructor.

NBK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945,
                         0.00600069867, 0.0200179144, 103.560653)
lens = SphericalLens(51.5mm, Inf, 3.6mm, BMO.inch, NBK7)   # plano-convex, like a stock f = 100 mm lens
translate3d!(lens, [0, 350mm, h2])
f = BMO.lensmakers_eq(51.5mm, Inf, NBK7(λ))               # ≈ 99.98 mm

# The first radius of curvature (51.5 mm) describes the curved front surface, which faces the collimated beam coming from the periscope; the second surface is flat (`Inf`).
#
# To model the coherent laser beam itself rather than a single ray, we use a [`GaussianBeamlet`](@ref) with a 0.5 mm waist radius. The `support` keyword pins the local reference frame of the beamlet to a fixed vector, which keeps the simulation reproducible instead of relying on an arbitrary default that could change from run to run.

w0 = 0.5mm
laser = GaussianBeamlet([0, 0, h1], [0, 1.0, 0], λ, w0; support=[1.0, 0, 0])
system = System([m1, m2, lens])
solve_system!(system, laser)

zs = range(500mm, 600mm, length=2001)      # optical path length from the laser
w, _, _, _ = gauss_parameters(laser, zs)
i = argmin(w)
println("waist: ", round(w[i] * 1e6, digits=1), " µm at z = ", round(zs[i] / mm, digits=1), " mm")

#-

fig3 = Figure(size=(600, 300))
ax3 = Axis(fig3[1,1], xlabel="optical path [mm]", ylabel="beam radius w [µm]")
lines!(ax3, zs ./ mm, w .* 1e6)
save("waist.png", fig3, px_per_unit=4); nothing #hide #md
fig3 #!md

# ![Beam radius as a function of optical path length, focused by the lens](waist.png)
#
# As a sanity check, we can compare this to the paraxial estimate ``\lambda f / (\pi w_{lens})``, using the beam radius at the lens (path length 450 mm). The beam has diverged considerably over that distance -- its Rayleigh range is only about 1.24 m -- so using the beam radius at the lens rather than the initial waist gives a more meaningful estimate here.

w_lens, _, _, _ = gauss_parameters(laser, 450mm)
println("paraxial waist estimate: ", round(λ * f / (π * w_lens) * 1e6, digits=1), " µm")

# The two values agree to within a few percent. Now we place a camera right at the waist location found above. The periscope path after `M2` runs along +y starting at `y = 100mm` (optical path 200 mm), so the camera's y-position follows from the path length at the waist.

camera = Detector(1mm)
translate3d!(camera, [0, 100mm + (zs[i] - 200mm), h2])
system = System([m1, m2, lens, camera])
laser = GaussianBeamlet([0, 0, h1], [0, 1.0, 0], λ, w0; support=[1.0, 0, 0])
empty!(camera)
solve_system!(system, laser)
x, z, I = intensity(camera; n=100)
println("power on camera: ", round(optical_power(camera) * 1e3, digits=3), " mW")

#-

fig4 = Figure(size=(500, 400))
ax4 = Axis(fig4[1,1], aspect=1, xlabel="x [µm]", ylabel="z [µm]")
hm = heatmap!(ax4, x .* 1e6, z .* 1e6, I)
Colorbar(fig4[1,2], hm)
save("camera.png", fig4, px_per_unit=4); nothing #hide #md
fig4 #!md

# ![Focused intensity distribution on the camera sensor](camera.png)
#
# Nearly all of the input power reaches the camera, and the focused spot matches the waist size computed above.

system_full = System([m1, m2, lens, camera])
figo = Figure(size=(600, 400))
axo = Axis3(figo[1,1], aspect=:data, azimuth=0.3π, elevation=0.15π)
hidedecorations!(axo)
render!(axo, system_full)
render!(axo, laser, color=:red)
save("laser_alignment.png", figo, px_per_unit=4); nothing #hide #md
figo #!md

# ## Next steps
#
# From here, you could continue with the [Michelson interferometer](@ref) tutorial to see how a similar beam is split and recombined, the [Miniature microscope](@ref) tutorial for a more complex multi-lens system, or browse the [Optical components](@ref) and [Visualization](@ref) sections for more details on the building blocks used here.
