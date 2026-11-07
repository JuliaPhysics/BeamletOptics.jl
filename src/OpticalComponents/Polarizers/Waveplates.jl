"""
    QuarterWavePlate(diameter; λ0 = nothing)

Creates a round, infinitesimally thin quarter-wave plate of the `diameter` in \\[m\\]: a
[`CoatingLayer`](@ref) with a [`Retarder`](@ref) of the retardance π/2. It faces the y-axis, and its fast
axis lies along the x-axis; rotate it about its normal with `yrotate3d!` to turn the fast axis, see
[`fast_axis`](@ref). With the design wavelength `λ0` in \\[m\\], it is a zero-order plate whose
retardance at the wavelength `λ` is `π/2 λ0/λ`, otherwise it is achromatic.

```julia
qwp = QuarterWavePlate(25.4e-3)
yrotate3d!(qwp, deg2rad(45))   # linear x into circular light
```
"""
QuarterWavePlate(diameter::Real; λ0::Nullable{Real} = nothing) = wave_plate(diameter, π / 2, λ0)

"""
    HalfWavePlate(diameter; λ0 = nothing)

Creates a round, infinitesimally thin half-wave plate: a [`CoatingLayer`](@ref) with a
[`Retarder`](@ref) of the retardance π. Turned by θ about its normal, it turns a linear polarization
along x by 2θ. Refer to [`QuarterWavePlate`](@ref) for its orientation and `λ0`.
"""
HalfWavePlate(diameter::Real; λ0::Nullable{Real} = nothing) = wave_plate(diameter, π, λ0)

# A round retarder film of the retardance Γ, fast axis along x
function wave_plate(diameter::Real, Γ::Real, λ0)
    diameter > 0 || throw(ArgumentError(lazy"the diameter of a wave plate must be positive, not $diameter"))
    return CoatingLayer(CircularFlatMesh(diameter / 2), Retarder(Γ, [1, 0, 0]; λ0))
end
