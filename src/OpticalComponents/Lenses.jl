"""
    AbstractRefractiveOptic <: AbstractObject

A generic type to represent an [`AbstractObject`](@ref) that refracts incoming rays.

# Implementation reqs.

Subtypes of `AbstractRefractiveOptic` should implement all supertype reqs. as well as:

## Fields

- `n`: a callable field which returns the [`RefractiveIndex`](@ref) for a wavelength `λ`

## Getters/setters

- `refractive_index`: gets the ref. index data of the optic

## Functions

- `interact3d`: the interaction logic should be akin to [`refraction3d`](@ref) for each surface crossing

# Additional information

The information provided below applies to the standard functional implementation of this type and may be overwritten
by specialized subtypes.

!!! info "Uniform ref. index"
    It is assumed that the optic consists of a single transparent material with a homogeneous refractive index `n`.
    It does not consider coated surfaces.

!!! info "Polarization ray tracing"
    Fresnel coefficients at the point of refraction are calculated via the [`fresnel_coefficients`](@ref) function with the
    refractive index data of the substrate and the previous medium.
"""
abstract type AbstractRefractiveOptic{T, S <: AbstractShape{T}, F} <: AbstractObject{T, S} end

refractive_index(object::AbstractRefractiveOptic) = object.n
refractive_index(object::AbstractRefractiveOptic{<:Any, <:Any, <:RefractiveIndex}, λ::Real)::Float64 = object.n(λ)

"""
    interact3d(AbstractSystem, AbstractRefractiveOptic, Beam, Ray)

Implements the refraction of a [`Ray`](@ref) at an optical surface. The "outside" ref. index is obtained from the `system` unless specified otherwise.
At the critical angle, total internal reflection occurs (see [`refraction3d`](@ref)).
"""
function interact3d(system::AbstractSystem,
        optic::AbstractRefractiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    # Check dir. of ray and surface normal
    normal = normal3d(intersection(ray))
    lambda = wavelength(ray)
    if isentering(ray)
        # Entering optic
        n1 = refractive_index(ray)
        n2 = refractive_index(optic, lambda)
        # Hint to test optic again
        hint = Hint(optic)
    else
        # Exiting optic
        n1 = refractive_index(optic, lambda)
        n2 = refractive_index(system, lambda)
        hint = nothing
        # Flip normal for refraction3d
        normal = -normal
    end
    # Calculate new dir. and pos.
    ndir, TIR = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    # In case of TIR, update hint and n2
    if TIR
        hint = Hint(optic)
        n2 = refractive_index(optic, lambda)
    end
    return BeamInteraction{T, R}(hint,
        Ray{T}(npos, ndir, nothing, wavelength(ray), n2))
end

"""
    interact3d(AbstractSystem, AbstractRefractiveOptic, Beam, PolarizedRay)

Implements the refraction of a [`PolarizedRay`](@ref) at an uncoated optical surface. The "outside" ref. index is obtained from the `system` unless specified otherwise.
Reflection and transmission values are calculated via the [`fresnel_coefficients`](@ref). Stray light is not tracked.
In the case of total internal reflection, only the reflected light is traced.
"""
function interact3d(system::AbstractSystem, optic::AbstractRefractiveOptic,
        ::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    lambda = wavelength(ray)
    normal = normal3d(intersection(ray))
    raypos = position(ray) + length(ray) * direction(ray)
    if isentering(ray)
        # Entering optic
        n1 = refractive_index(ray)
        n2 = refractive_index(optic, lambda)
        # Hint to test optic again
        hint = Hint(optic)
    else
        # Exiting optic
        n1 = refractive_index(optic, lambda)
        n2 = refractive_index(system, lambda)
        hint = nothing
        # Flip normal for refraction3d
        normal = -normal
    end
    # Calculate (and correct into 1. quadrant) the angle of incidence
    θi = angle3d(direction(ray), -normal)
    # Get Fresnel coefficients
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    # Optical interaction
    if is_internally_reflected(rp, rs)
        # Update hint and outgoing ref. index
        hint = Hint(optic)
        n2 = refractive_index(optic, lambda)
        # Calculate reflection
        new_dir = reflection3d(direction(ray), normal)
        J = [-rs 0 0; 0 rp 0; 0 0 1]
    else
        # Calculate refraction
        new_dir, ~ = refraction3d(direction(ray), normal, n1, n2)
        J = [ts 0 0; 0 tp 0; 0 0 1]
    end
    # Calculate new polarization
    E0 = _calculate_global_E0(direction(ray), new_dir, J, polarization(ray))
    return BeamInteraction{T, R}(
        hint, PolarizedRay{T}(raypos, new_dir, nothing, wavelength(ray), n2, E0))
end

"""
    Lens{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, S, N}

Represents an uncoated `Lens` with a homogeneous [`RefractiveIndex`](@ref) `n = n(λ)`.
Refer to the [`Lens`](@ref) and [`SphericalLens`](@ref) constructors for more information on how to generate lenses.

# Fields

- `shape`: geometry of the lens, refer to [`AbstractShape`](@ref) for more information
- `n`: [`RefractiveIndex`](@ref) function that returns n(λ)

# Additional information

!!! info "Refractive index"
    The chromatic dispersion of the lens is represented by a λ-dependent function for `n`
    and must be provided by the user. For testing purposes, an anonymous function, e.g. λ -> 1.5
    can be passed such that the lens has the same refractive index for all wavelengths.
"""
struct Lens{T, S <: AbstractShape{T}, N <: RefractiveIndex} <:
       AbstractRefractiveOptic{T, S, N}
    shape::S
    n::N
    function Lens(
            shape::S, n::N) where {T <: Real, S <: AbstractShape{T}, N <: RefractiveIndex}
        test_refractive_index_function(n)
        return new{T, S, N}(shape, n)
    end
end

thickness(l::Lens) = thickness(shape(l))

"""
     Lens(front_surface::AbstractRotationallySymmetricSurface, back_surface::AbstractRotationallySymmetricSurface, center_thickness::Real, n::RefractiveIndex)

Constructs a new [`Lens`](@ref) object using the surface specifications `front_surface` and
`back_surface` and the `center_thickness`. These inputs are used to construct a [`UnionSDF`](@ref)
that consists of the appropriate sub-SDFs to represent the shape of the lens.

The material properties are supplied via the `n` parameter.

# Additional information

!!! info "Radius of curvature (ROC) sign definition"
    The ROC is defined to be positive if the center is to the right of the surface. Otherwise it is negative.

!!! warning "Meniscus"
    If your specification results in a meniscus lens, only spherical meniscus lenses are supported at the moment.
"""
function Lens(
        front_surface::AbstractRotationallySymmetricSurface,
        back_surface::AbstractRotationallySymmetricSurface,
        center_thickness::Real,
        n::RefractiveIndex)
    # Define effective (optical and mechanical) diameters:
    d_mid = min(diameter(front_surface), diameter(back_surface))
    md_mid = max(mechanical_diameter(front_surface), mechanical_diameter(back_surface))

    # Initialize remaining cylindrical section length.
    l0 = center_thickness

    # Front Surface
    front = sdf(front_surface, ForwardOrientation())
    l0 -= isnothing(front) ? zero(l0) : thickness(front)

    # Back Surface
    back = sdf(back_surface, BackwardOrientation())
    l0 -= isnothing(back) ? zero(l0) : thickness(back)

    # Use MeniscusLensSDF if cylinder length is non-positive
    if l0 ≤ 0
        if sign(radius(front_surface)) == sign(radius(back_surface))
            shape = meniscus_lens_sdf(
                front_surface, front, back_surface, back, center_thickness)
            # add an outer ring if necessary
            if md_mid > d_mid
                _thickness = thickness(shape.cylinder)
                pos = position(shape.cylinder)
                ring = RingSDF(d / 2, (md_mid - d_mid) / 2, _thickness)
                translate3d!(ring, [0, pos[2] + _thickness / 2, 0])
                shape += ring
            end
        else
            throw(ArgumentError("Lens parameters lead to cylinder section length of ≤ 0, use ThinLens instead."))
        end
    else
        # Construct central plano surface and add front/back surfaces
        mid = PlanoSurfaceSDF(l0, d_mid)
        if front !== nothing
            translate3d!(mid, [0, thickness(front), 0])
            mid += front
        end
        if back !== nothing
            translate3d!(back, [0, thickness(mid) + thickness(back), 0])
            mid += back
        end
        shape = mid

        # Add mechanical ring if md_mid > d_mid
        if md_mid > d_mid
            ring_thickness = thickness(mid)
            ring_center = position(mid)[2] + ring_thickness / 2
            if front !== nothing
                s = edge_sag(front_surface, front)
                ring_thickness -= s
                ring_center += s / 2
            end
            if back !== nothing
                s = edge_sag(back_surface, back)
                ring_thickness += s
                ring_center += s / 2
            end
            ring = RingSDF(d_mid / 2, (md_mid - d_mid) / 2, ring_thickness)
            translate3d!(ring, [0, ring_center, 0])
            shape += ring
        elseif md_mid < d_mid
            @warn "Mechanical diameter is less than clear aperture; parameter md has been ignored."
        end
    end

    return Lens(shape, n)
end

function Lens(
        front_surface::Union{Nothing, AbstractCylindricSurface},
        back_surface::Union{Nothing, AbstractCylindricSurface},
        center_thickness::Real,
        n::RefractiveIndex;
        circular=false,
        circular_radius=Inf)
    # Initialize remaining box section length.
    l0 = center_thickness

    # Front Surface
    front = isnothing(front_surface) ? nothing : sdf(front_surface, ForwardOrientation())
    l0 -= isnothing(front) ? zero(l0) : thickness(front)

    # Back Surface
    back = isnothing(back_surface) ? nothing : sdf(back_surface, BackwardOrientation())
    l0 -= isnothing(back) ? zero(l0) : thickness(back)

    # validate
    if front_surface !== nothing && back_surface !== nothing
        height(front_surface) != height(back_surface) && throw(ArgumentError("height of front and back surface have to match for cylindric lenses"))

        d_mid = min(diameter(front_surface), diameter(back_surface))
        md_mid = max(mechanical_diameter(front_surface), mechanical_diameter(back_surface))
        h = height(front_surface)
    elseif front_surface !== nothing
        d_mid = diameter(front_surface)
        md_mid = mechanical_diameter(front_surface)
        h = height(front_surface)
    elseif back_surface !== nothing
        d_mid = diameter(back_surface)
        md_mid = mechanical_diameter(back_surface)
        h = height(back_surface)
    else
        throw(ArgumentError("Both surfaces cannot be nothing"))
    end

    # Check if the center section is positive
    if l0 ≤ 0
        throw(ArgumentError("Lens parameters lead to a box section length of ≤ 0"))
    else
        # Construct central box and add front/back surfaces
        mid = BoxSDF(h, l0, d_mid)
        translate3d!(mid, [0, l0/2, 0])
        if front !== nothing
            translate3d!(mid, [0, thickness(front), 0])
            mid += front
        end
        if back !== nothing
            translate3d!(back, [0, thickness(mid) + thickness(back), 0])
            mid += back
        end
        shape = mid

        # Add mechanical ring if md_mid > d_mid
        if md_mid > d_mid
            throw(ArgumentError("Mechanical diameters not yet implemented for cylindric lenses"))
            # ring_thickness = thickness(mid)
            # ring_center = position(mid)[2] + ring_thickness / 2
            # if front !== nothing
            #     s = edge_sag(front_surface, front)
            #     ring_thickness -= s
            #     ring_center += s / 2
            # end
            # if back !== nothing
            #     s = edge_sag(back_surface, back)
            #     ring_thickness += s
            #     ring_center += s / 2
            # end
            # ring = RingSDF(d_mid / 2, (md_mid - d_mid) / 2, ring_thickness)
            # translate3d!(ring, [0, ring_center, 0])
            # shape += ring
        elseif md_mid < d_mid
            @warn "Mechanical diameter is less than clear aperture; parameter md has been ignored."
        end
    end

    return Lens(shape, n)
end
