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
abstract type AbstractRefractiveOptic{T, F} <: AbstractObject{T} end

refractive_index(object::AbstractRefractiveOptic) = object.n
refractive_index(object::AbstractRefractiveOptic{<:Any, <:RefractiveIndex}, λ::Real)::Float64 = object.n(λ)

medium(optic::AbstractRefractiveOptic) = IsotropicMedium(refractive_index(optic))

"""
    interact3d(AbstractSystem, AbstractRefractiveOptic, Beam, AbstractRay)

Implements the interaction of a [`Ray`](@ref) or a [`PolarizedRay`](@ref) at a surface of a refracting
optic: the [`transition`](@ref) at the hit gives the media on both sides, the optic, the system
medium, or another optic that touches the surface, e.g. the other lens of a cemented doublet. The
[`surface_model`](@ref) of the boundary acts between them. For an uncoated surface, the ray is
refracted, or totally reflected beyond the critical angle (see [`cross_boundary`](@ref)); a
[`PolarizedRay`](@ref) is transmitted with the [`fresnel_coefficients`](@ref), stray light is not tracked.
"""
function interact3d(system::AbstractSystem,
        ::AbstractRefractiveOptic,
        beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    t = transition(system, ray)
    return interact3d(system, surface_model(t), t, beam, ray)
end

"""
    Lens{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, N}

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
struct Lens{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, N}
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
    _signed_edge_sag(surface, sd)

Returns the edge sagitta of the `surface` with sign, such that the edge of a front surface lies at
`y = s` and the edge of a back surface at `y = center_thickness + s`. [`edge_sag`](@ref) returns
the unsigned sagitta for spherical and cylindrical surfaces, hence the sign of the radius is
applied; the sagitta of aspherical and acylindrical surfaces is already signed.
"""
_signed_edge_sag(s::AbstractSurface, sd) = edge_sag(s, sd)
_signed_edge_sag(s::SphericalSurface, sd) = sign(radius(s)) * edge_sag(s, sd)
_signed_edge_sag(s::CylindricalSurface, sd) = sign(radius(s)) * edge_sag(s, sd)

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
    A meniscus lens, whose center thickness does not exceed the sagitta of its convex surface, must consist of
    spherical surfaces. An aspherical surface throws an `ArgumentError` in this case.
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
                ring = RingSDF(d_mid / 2, (md_mid - d_mid) / 2, _thickness)
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

        d_front = diameter(front_surface)
        d_back = diameter(back_surface)
        d_min = min(d_front, d_back)
        d_max = max(d_front, d_back)

        if md_mid < d_min
            @warn "Mechanical diameter is less than clear aperture; parameter md has been ignored."
        else
            # add leveling ring if there is a step between front and back clear apertures
            if d_front != d_back
                if d_back > d_front
                    # Step exists on the front side: level front to match back.
                    leveling_thickness = l0
                    if front !== nothing
                        s_front = _signed_edge_sag(front_surface, front)
                        if s_front < 0
                            # edge curves towards negative so it is a concave type shape,
                            # which has to be covered by the ring
                            leveling_thickness += abs(s_front) + thickness(front)
                        end
                    end

                    ring = RingSDF(d_front / 2, (d_back - d_front) / 2, leveling_thickness)
                    translate3d!(ring, [0, _signed_edge_sag(front_surface, front) + leveling_thickness / 2, 0])
                    shape += ring
                else  # d_front > d_back
                    # Step exists on the back side: level back to match front.
                    leveling_thickness = l0
                    if back !== nothing
                        s_back = _signed_edge_sag(back_surface, back)
                        if (s_back - thickness(back)) > 0
                            # edge curves towards positive so it is a concave type shape,
                            # which has to be covered by the ring
                            leveling_thickness += abs(s_back) + thickness(back)
                        end
                    end
                    leveling_center = position(mid)[2] +
                                      (l0/2 + (back !== nothing ? _signed_edge_sag(back_surface, back) : 0))
                    if back !== nothing
                        leveling_center += s_back / 2
                    end
                    ring = RingSDF(d_back / 2, (d_front - d_back) / 2, leveling_thickness)
                    translate3d!(ring, [0, thickness(front) + leveling_thickness/2, 0])
                    shape += ring
                end
            end

            # add outer ring if the mechanical diameter exceeds the larger clear aperture.
            if md_mid > d_max
                outer_thickness = thickness(mid)
                outer_center = position(mid)[2] + outer_thickness / 2
                if front !== nothing
                    s_front = _signed_edge_sag(front_surface, front)
                    outer_thickness -= s_front
                    outer_center += s_front / 2
                end
                if back !== nothing
                    s_back = _signed_edge_sag(back_surface, back)
                    outer_thickness += s_back
                    outer_center += s_back / 2
                end
                ring = RingSDF(d_max / 2, (md_mid - d_max) / 2, outer_thickness)
                translate3d!(ring, [0, outer_center, 0])
                shape += ring
            end
        end
    end

    return Lens(shape, n)
end

function Lens(front_surface::AbstractRotationallySymmetricSurface,
        center_thickness::Real, n::RefractiveIndex)
    Lens(
        front_surface,
        CircularFlatSurface(diameter(front_surface)),
        center_thickness,
        n
    )
end

function Lens(front_surface::CircularFlatSurface, back_surface::CircularFlatSurface,
        center_thickness::Real, n::RefractiveIndex)
    d_mid = min(diameter(front_surface), diameter(back_surface))

    return Lens(
        PlanoSurfaceSDF(center_thickness, d_mid),
        n
    )
end

"""
     Lens(front_surface::AbstractCylindricalSurface, back_surface::AbstractCylindricalSurface, center_thickness::Real, n::RefractiveIndex)

Constructs a new [`Lens`](@ref) object using the cylindric surface specifications `front_surface` and
`back_surface` and the `center_thickness`. These inputs are used to construct a [`UnionSDF`](@ref)
that consists of the appropriate sub-SDFs to represent the shape of the lens.

This method of `Lens` is specific for cylindric lenses and has some limitations:
    - The cylinder height of both surfaces has to be identical
    - No mixture with non-cylindric surfaces is supported at the moment

The material properties are supplied via the `n` parameter.

# Additional information

!!! info "Radius of curvature (ROC) sign definition"
    The ROC is defined to be positive if the center is to the right of the surface. Otherwise it is negative.
"""
function Lens(
        front_surface::AbstractCylindricalSurface,
        back_surface::AbstractCylindricalSurface,
        center_thickness::Real,
        n::RefractiveIndex)
    # Initialize remaining box section length.
    l0 = center_thickness

    # Front Surface
    front = sdf(front_surface, ForwardOrientation())
    l0 -= isnothing(front) ? zero(l0) : thickness(front)

    # Back Surface
    back = sdf(back_surface, BackwardOrientation())
    l0 -= isnothing(back) ? zero(l0) : thickness(back)

    # validate
    d_mid, md_mid, h = cylindric_lens_outer_parameters(front_surface, back_surface)

    # Check if the center section is positive
    if l0 ≤ 0
        throw(ArgumentError("Lens parameters lead to a box section length of ≤ 0"))
    else
        # Construct central box and add front/back surfaces
        mid = BoxSDF(h, l0, d_mid)
        translate3d!(mid, [0, l0 / 2, 0])
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
                s = _signed_edge_sag(front_surface, front)
                ring_thickness -= s
                ring_center += s / 2
            end
            if back !== nothing
                s = _signed_edge_sag(back_surface, back)
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

function Lens(front_surface::AbstractCylindricalSurface,
        center_thickness::Real, n::RefractiveIndex)
    Lens(
        front_surface,
        RectangularFlatSurface(diameter(front_surface)),
        center_thickness,
        n
    )
end

function cylindric_lens_outer_parameters(
        f::AbstractCylindricalSurface, b::AbstractCylindricalSurface)
    height(f) != height(b) &&
        throw(ArgumentError("height of front and back surface have to match for cylindric lenses"))

    d_mid = min(diameter(f), diameter(b))
    md_mid = max(mechanical_diameter(f), mechanical_diameter(b))
    h = height(f)

    return d_mid, md_mid, h
end

function cylindric_lens_outer_parameters(
        f::AbstractCylindricalSurface, ::RectangularFlatSurface)
    d_mid = diameter(f)
    md_mid = mechanical_diameter(f)
    h = height(f)

    return d_mid, md_mid, h
end

function cylindric_lens_outer_parameters(
        f::RectangularFlatSurface, b::AbstractCylindricalSurface)
    cylindric_lens_outer_parameters(b, f)
end

function cylindric_lens_outer_parameters(
        f::RectangularFlatSurface, b::RectangularFlatSurface)
    d_mid = diameter(f)
    md_mid = mechanical_diameter(f)

    return d_mid, md_mid, d_mid
end

function Lens(front_surface::RectangularFlatSurface, back_surface::RectangularFlatSurface,
        center_thickness::Real, n::RefractiveIndex)
    d_mid = min(diameter(front_surface), diameter(back_surface))
    mid = BoxSDF(d_mid, center_thickness, d_mid)
    translate3d!(mid, [0, center_thickness / 2, 0])

    return Lens(
        mid,
        n
    )
end
