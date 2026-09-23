using GLMakie, BeamletOptics

const BMO = BeamletOptics
const cm = 1e-2
const mm = 1e-3

## Michelson interferometer, see the Michelson tutorial
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

## Figure with the 3D view, the detector intensity and the optical power
fig = Figure(size = (1200, 700))
ax = LScene(fig[1:2, 1]; show_axis = false)
heat_ax = Axis(fig[1, 2]; title = "Detector intensity", xlabel = "x [mm]", ylabel = "y [mm]", aspect = 1)
power_ax = Axis(fig[2, 2]; title = "Optical power", xlabel = "Update", ylabel = "P [mW]")
status = Label(fig[3, 1:2], "Click on a component to select it, press h to show the controls"; tellwidth = false)
colsize!(fig.layout, 1, Relative(0.6))

system_handle = live_render!(ax, system)
beam_handle = live_render!(ax, beam)

# Evaluate the detector on a fixed grid, the color range is set by the aligned interferometer
n = 100
detector_intensity() = intensity(pd; n, x_min = -pd_size / 2, x_max = pd_size / 2, z_min = -pd_size / 2, z_max = pd_size / 2)
x, y, I = detector_intensity()
dA = step(x) * step(y)
I_obs = Observable(I)
heatmap!(heat_ax, x / mm, y / mm, I_obs; colorrange = (0, maximum(I)))

power = Observable([Point2f(1, 1e3 * sum(I) * dA)])
lines!(power_ax, power; color = :red)

## Solve the system and update all plots after a component has been moved
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

# Open the interactive window when used from the REPL or run as a script
if isinteractive()
    display(fig)
elseif abspath(PROGRAM_FILE) == @__FILE__
    wait(display(fig))
end
