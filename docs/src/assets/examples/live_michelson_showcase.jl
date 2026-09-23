include(joinpath(@__DIR__, "live_michelson.jl"))

GLMakie.activate!(; visible = false)

## Tilt mirror 1 by 1 mrad, which generates fringes on the detector
controls.selected[] = m1
for _ in 1:10
    zrotate3d!(m1, 1e-4)
    on_change(m1)
end
save("live_michelson_fringes.png", fig; px_per_unit = 2)

## Undo the tilt and move mirror 2 by λ/2, i.e. one period of the optical power
zrotate3d!(m1, -1e-3)
controls.selected[] = m2
for _ in 1:40
    translate3d!(m2, [0, λ / 80, 0])
    on_change(m2)
end
save("live_michelson_power.png", fig; px_per_unit = 2)
