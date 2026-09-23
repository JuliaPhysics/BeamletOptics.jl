```@setup live_michelson
dir = joinpath(@__DIR__, "..", "assets", "examples")

Main.DocUtils.conditional_include(joinpath(dir, "live_michelson_showcase.jl"))
```

# Interactive Michelson interferometer

This example shows how to build an interactive application with the [live rendering](@ref "Live rendering") functions of this package. The Michelson interferometer of the [Michelson interferometer](@ref) tutorial is rendered into a `GLMakie` window, in which the components can be moved and rotated with the mouse and keyboard. After each change, the system is solved again and the beam path, the fringe pattern on the detector and the optical power are updated live.

![Interactive Michelson interferometer](live_michelson_fringes.png)

The full script can be found [here](https://github.com/JuliaPhysics/BeamletOptics.jl/blob/master/docs/src/assets/examples/live_michelson.jl). Running it from the REPL opens the interactive window.

## Setting up the system

The system is identical to the tutorial, but without the optomechanical parts. A [`GaussianBeamlet`](@ref) is used as the HeNe laser source.

```julia
using GLMakie, BeamletOptics

const BMO = BeamletOptics
const cm = 1e-2
const mm = 1e-3

λ = 632.8e-9
w0 = 0.65mm / 2
M2 = 0.7e-3 / BMO.divergence_angle(λ, w0, 1)
beam = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; M2)

NBK7 = DiscreteRefractiveIndex([λ], [1.51509])

rpm = RightAnglePrismMirror(25mm, 25mm)
zrotate3d!(rpm, deg2rad(45))
translate3d!(rpm, [0, 23.5cm, 0])

cbs = CubeBeamsplitter(BMO.inch, NBK7)
zrotate3d!(cbs, deg2rad(-90))
translate_to3d!(cbs, [18.81cm, 23.5cm, 0])

m1 = RoundPlanoMirror(BMO.inch, 5mm)
zrotate3d!(m1, deg2rad(-90))
translate3d!(m1, [42.715cm, 23.5cm, 0])

m2 = RoundPlanoMirror(BMO.inch, 5mm)
translate3d!(m2, [18.81cm, 37.405cm, 0])

pd_size = 5mm
pd = Detector(pd_size)
translate_to3d!(pd, [18.81cm, 9.595cm, 0])

system = System([rpm, cbs, m1, m2, pd])
solve_system!(system, beam)
```

## Rendering the system

The figure consists of a `LScene` for the system and two axes for the detector intensity and the optical power. Instead of [`render!`](@ref), the system and the beam are rendered with [`live_render!`](@ref), which returns handles that can be updated later on.

```julia
fig = Figure(size = (1200, 700))
ax = LScene(fig[1:2, 1]; show_axis = false)
heat_ax = Axis(fig[1, 2]; title = "Detector intensity", xlabel = "x [mm]", ylabel = "y [mm]", aspect = 1)
power_ax = Axis(fig[2, 2]; title = "Optical power", xlabel = "Update", ylabel = "P [mW]")
status = Label(fig[3, 1:2], "Click on a component to select it, press h to show the controls"; tellwidth = false)
colsize!(fig.layout, 1, Relative(0.6))

system_handle = live_render!(ax, system)
beam_handle = live_render!(ax, beam)
```

The detector intensity is evaluated on a fixed grid and stored in an `Observable`, such that the heatmap follows its changes. The color range is set by the aligned interferometer, which makes changes of the optical power visible.

```julia
n = 100
detector_intensity() = intensity(pd; n, x_min = -pd_size / 2, x_max = pd_size / 2, z_min = -pd_size / 2, z_max = pd_size / 2)
x, y, I = detector_intensity()
dA = step(x) * step(y)
I_obs = Observable(I)
heatmap!(heat_ax, x / mm, y / mm, I_obs; colorrange = (0, maximum(I)))

power = Observable([Point2f(1, 1e3 * sum(I) * dA)])
lines!(power_ax, power; color = :red)
```

## Interaction

The `on_change` function is called after a component has been moved. It solves the system again, updates the beam via [`update_render!`](@ref) and evaluates the detector. The rendering of the moved component itself is updated automatically. Note that the detector has to be reset via `empty!` before solving the system.

```julia
function on_change(obj)
    empty!(pd)
    solve_system!(system, beam)
    update_render!(beam_handle)
    # No light hits the detector if e.g. the beamsplitter is moved out of the beam
    I_obs[] = isnothing(BMO.hits(pd)) ? zero(I_obs[]) : detector_intensity()[3]
    P = sum(I_obs[]) * dA
    push!(power[], Point2f(last(power[])[1] + 1, 1e3 * P))
    length(power[]) > 300 && popfirst!(power[])
    notify(power)
    autolimits!(power_ax)
    status.text[] = "$(nameof(typeof(obj))) at $(round.(BMO.position(obj) / mm, digits = 6)) mm, P = $(round(1e3 * P, digits = 4)) mW"
    return nothing
end

controls = kinematic_controls!(ax, system_handle; on_change)
display(fig)
```

[`kinematic_controls!`](@ref) enables the following controls:

| Input                              | Effect                                                     |
|:-----------------------------------|:-----------------------------------------------------------|
| Left-drag on a component           | Move it in the horizontal plane                            |
| Shift + left-drag on a component   | Rotate it around the vertical axis                         |
| `↑` / `↓`                          | Move it along its normal by 10 nm                          |
| `←` / `→`                          | Rotate it around the vertical axis by 10 µrad              |
| `Page Up` / `Page Down`            | Tilt it around its local x-axis by 10 µrad                 |
| Shift + key                        | Ten times the step size                                    |
| `r`                                | Reset the component to its initial pose                    |
| `Esc` / click on empty space       | Deselect the component                                     |
| `h`                                | Show or hide an overlay of all controls                    |

The selected component is marked by a box and three arrows, which show the direction of `↑` (green), the tilt axis of `Page Up` (red) and the rotation axis of `←` (blue). 

Since the interferometer is sensitive to changes in the order of the wavelength, the keyboard controls are best suited for alignment. Rotating the mirror `m1` by 1 mrad generates the fringes shown above. Moving the mirror `m2` by ``\lambda/2`` changes the optical path length by ``\lambda``, which corresponds to one period of the optical power:

![Optical power](live_michelson_power.png)
