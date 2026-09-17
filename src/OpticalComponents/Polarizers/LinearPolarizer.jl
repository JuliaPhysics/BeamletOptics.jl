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

"""
    _is_film_interface(lp::LinearPolarizer, from, ray::AbstractRay)::Bool

Tests whether `ray`, having just crossed a prism (`from` is `lp.front` or `lp.back`), is exiting through
the flat inner face that is cemented to the polarizing film (as opposed to exiting the outer face or the
cylindrical edge close to the film).
"""
function _is_film_interface(lp::LinearPolarizer, from, ray::AbstractRay)
    isentering(ray) && return false
    y = orientation(lp)[:, 2]
    nrm = normal3d(intersection(ray))
    abs(dot(nrm, y)) ≥ 0.5 || return false
    c = position(lp)
    p = position(ray) + length(ray) * direction(ray)
    if from === lp.front
        return dot(p - c, y) > -thickness(lp.front) / 2
    else
        return dot(p - c, y) < thickness(lp.back) / 2
    end
end

"""
    _film_interface(system, lp::LinearPolarizer, from, to, beam, ray::Ray)

Handles the cemented `from → to` glass interface for an unpolarized [`Ray`](@ref): the film has no effect,
only glass-to-glass refraction (with total internal reflection handling) occurs.
"""
function _film_interface(system::AbstractSystem, lp::LinearPolarizer, from, to,
        beam::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    λ = wavelength(ray)
    n1 = refractive_index(from, λ)
    n2 = refractive_index(to, λ)
    normal = -normal3d(intersection(ray))
    ndir, TIR = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    if TIR
        hint = Hint(lp, shape(from))
        n = n1
    else
        hint = Hint(lp, shape(to))
        n = n2
    end
    return BeamInteraction{T, R}(hint, Ray{T}(npos, ndir, nothing, λ, n))
end

"""
    _film_interface(system, lp::LinearPolarizer, from, to, beam, ray::PolarizedRay)

Handles the cemented `from → to` glass interface for a [`PolarizedRay`](@ref). The polarizing film acts on
the incident field (still expressed in `from`) before the ray crosses into `to`. In case of total internal
reflection at the cemented interface, the ray is reflected back into `from` and the film is not applied.
"""
function _film_interface(system::AbstractSystem, lp::LinearPolarizer, from, to,
        beam::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    λ = wavelength(ray)
    n1 = refractive_index(from, λ)
    n2 = refractive_index(to, λ)
    normal = -normal3d(intersection(ray))
    npos = position(ray) + length(ray) * direction(ray)
    in_dir = direction(ray)
    θi = angle3d(in_dir, -normal)
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    if is_internally_reflected(rp, rs)
        new_dir = reflection3d(in_dir, normal)
        J = SPBasis(-rs, 0, 0, rp)
        E = polarization(ray)
        hint = Hint(lp, shape(from))
        n = n1
    else
        E = _calculate_global_E0(lp.filter, ray, in_dir, lp.filter.JMat)
        if norm(E) ≤ lp.filter.cutoff
            return nothing
        end
        new_dir, ~ = refraction3d(in_dir, normal, n1, n2)
        J = SPBasis(ts, 0, 0, tp)
        hint = Hint(lp, shape(to))
        n = n2
    end
    P = _calculate_global_E0(in_dir, new_dir, normal3d(intersection(ray)), J)
    E0 = P * E
    return BeamInteraction{T, R}(hint, PolarizedRay{T}(npos, new_dir, nothing, λ, n, E0))
end

"""
    interact3d(system::AbstractSystem, lp::LinearPolarizer, beam::Beam, ray::AbstractRay)

Dispatches the optical interaction between a [`Beam`](@ref)/[`AbstractRay`](@ref) and a [`LinearPolarizer`](@ref)
depending on which sub-shape (`front`, `filter` or `back`) was intersected.

# Additional information

!!! info "Cemented interface"
    The `front` and `back` substrates sit flush against the film. Refraction at the outer, uncoated surfaces is
    handled by the standard [`AbstractRefractiveOptic`](@ref) interaction logic. When a ray reaches the flat inner
    face cemented to the film, the film's Jones matrix is applied to the field before it is refracted (or, in case
    of total internal reflection, reflected without crossing the film) directly from the current substrate's
    refractive index into the other substrate's.
"""
function interact3d(system::AbstractSystem, lp::LinearPolarizer,
        beam::Beam{T, R}, ray::R) where {T <: Real, R <: AbstractRay{T}}
    hit = shape(intersection(ray))
    if hit === shape(lp.front) || hit === shape(lp.back)
        from, to = hit === shape(lp.front) ? (lp.front, lp.back) : (lp.back, lp.front)
        _is_film_interface(lp, from, ray) && return _film_interface(system, lp, from, to, beam, ray)
        i = interact3d(system, from, beam, ray)
        # Re-route the substrate's own hint through the polarizer, otherwise the next
        # hit is dispatched to the bare prism and the film interface is skipped
        if !isnothing(i) && !isnothing(hint(i))
            hint!(i, Hint(lp, shape(hint(i))))
        end
        return i
    elseif hit === shape(lp.filter)
        return interact3d(system, lp.filter, beam, ray)
    end
    error("LinearPolarizer: intersected shape is not part of this polarizer")
end
