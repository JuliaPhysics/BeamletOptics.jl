"""
    LinearPolarizer{T, N <: RefractiveIndex} <: AbstractObject{T}

Represents a real, round linear polarizing film laminated between two uncoated glass plates
(e.g. Thorlabs LPNIRE100-B), modeled as a thin [`PolarizationFilter`](@ref) cemented between two
[`Prism`](@ref) substrates that sit flush against the film. Refer to [`RoundLinearPolarizer`](@ref)
for the constructor.

# Fields

- `filter`: the thin [`PolarizationFilter`](@ref) that models the polarizing film
- `front`: front glass substrate [`Prism`](@ref)
- `back`: back glass substrate [`Prism`](@ref)

# Additional information

!!! info "Kinematic center"
    The center of kinematics of this component lies at the center of the film.

!!! info "Uncoated surfaces"
    Since the glass substrates are uncoated, a [`PolarizedRay`](@ref) loses approx. 4 % of its intensity per
    outer surface at n ≈ 1.5 due to Fresnel reflection, unlike the anti-reflection coated real-world part.

!!! info "Crossed orientation"
    When the film orientation is crossed (transmission axis perpendicular to the incident polarization),
    the beam is terminated at the cemented film interface once the norm of the transmitted electric field
    falls to or below the `cutoff_strength` of the underlying [`PolarizationFilter`](@ref).
"""
struct LinearPolarizer{T, N <: RefractiveIndex} <: AbstractObject{T}
    filter::PolarizationFilter{T, Mesh{T}}
    front::Prism{T, PlanoSurfaceSDF{T}, N}
    back::Prism{T, PlanoSurfaceSDF{T}, N}
end

shape_trait_of(::LinearPolarizer) = MultiShape()
# filter first: MultiShape uses the first element as kinematic center (pivot at film center)
shape(lp::LinearPolarizer) = (lp.filter, lp.front, lp.back)
refractive_index(lp::LinearPolarizer, λ::Real) = refractive_index(lp.front, λ)
thickness(lp::LinearPolarizer) = dot(position(lp.back) - position(lp.front), orientation(lp)[:, 2]) + thickness(lp.back)

"""
    transmission_axis(lp::LinearPolarizer)

Returns the unit vector (in global coordinates) along which [`LinearPolarizer`](@ref) `lp` transmits
polarization. Refer to [`transmission_axis(::PolarizationFilter)`](@ref) for details.
"""
transmission_axis(lp::LinearPolarizer) = transmission_axis(lp.filter)

"""
    RoundLinearPolarizer(diameter, front_thickness, back_thickness, n; cutoff_strength=eps())

Creates a [`LinearPolarizer`](@ref): a round polarizing film cemented between two round glass plates
of the given `diameter`, flush against the film on both sides.

# Inputs

- `diameter`: outer diameter of the film and substrates in [m]
- `front_thickness`: thickness of the front glass substrate in [m]
- `back_thickness`: thickness of the back glass substrate in [m]
- `n`: [`RefractiveIndex`](@ref) of both glass substrates

# Keywords

- `cutoff_strength`: passed to [`RoundPolarizationFilter`](@ref)

# Additional information

The component is centered on the film (local origin); the front substrate extends from
`y = -front_thickness` to `y = 0`, the back substrate from `y = 0` to `y = back_thickness`,
along local +y. The film transmits along local x and blocks local z.
"""
function RoundLinearPolarizer(diameter::Real, front_thickness::Real, back_thickness::Real,
        n::RefractiveIndex; cutoff_strength = eps())
    front_thickness > 0 || throw(ArgumentError("front_thickness must be positive"))
    back_thickness > 0 || throw(ArgumentError("back_thickness must be positive"))

    diameter, front_thickness, back_thickness = float.((diameter, front_thickness, back_thickness))

    filter = RoundPolarizationFilter(diameter; cutoff_strength)
    front = Prism(PlanoSurfaceSDF(front_thickness, diameter), n)
    translate3d!(front, [0, -front_thickness, 0])
    back = Prism(PlanoSurfaceSDF(back_thickness, diameter), n)
    return LinearPolarizer(filter, front, back)
end
