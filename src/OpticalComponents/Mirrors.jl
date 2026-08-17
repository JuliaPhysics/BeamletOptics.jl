"""
    AbstractReflectiveOptic <: AbstractObject

A generic type to represent an [`AbstractObject`] which reflects incoming rays.

# Implementation reqs.

Subtypes of `AbstractReflectiveOptic` should implement all supertype reqs. as well as:

## Fields

- no specific fields required

## Getters/setters

- none required

## Functions

- `interact3d`:  the interaction logic should be akin to [`reflection3d`](@ref) for each surface crossing

# Additional information

The information provided below applies to the standard functional implementation of this type and may be overwritten
by specialized subtypes.

!!! info "Polarization ray tracing"
    Fresnel coefficients during reflection are set such that no reflection losses occur (i.e. `|rₚ| = |rₛ| = 1`).
"""
abstract type AbstractReflectiveOptic{T} <: AbstractObject{T} end

@inline function resolve_coincident_boundary(
        exiting_int, entering_int, ::AbstractObject, ::AbstractReflectiveOptic)
    entering_int.coincident_object = object(exiting_int)
    return entering_int
end

@inline function resolve_coincident_boundary(
        exiting_int, entering_int, ::AbstractReflectiveOptic, ::AbstractObject)
    exiting_int.coincident_object = object(entering_int)
    return exiting_int
end

@inline function resolve_coincident_boundary(
        exiting_int, entering_int, ::AbstractReflectiveOptic, ::AbstractReflectiveOptic)
    exiting_int.coincident_object = object(entering_int)
    return exiting_int
end

function interact3d(system::AbstractSystem,
        optic::AbstractReflectiveOptic,
        beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    coated_obj, coating = resolve_coated_boundary(system, optic, ray)
    target_obj = isnothing(coated_obj) ? optic : coated_obj
    behavior = coating_behavior(coating, ray)
    if behavior isa Absorptive
        return nothing
    elseif !(coating isa Uncoated) && (behavior isa Transmissive || behavior isa Splitting)
        return interact_refractive_boundary(system, target_obj, coating, beam, ray)
    end
    return interact_reflective_boundary(system, target_obj, coating, beam, ray)
end

"""
    Mirror{S <: AbstractShape} <: AbstractReflectiveOptic

Concrete implementation of a perfect mirror (R = 1) with arbitrary shape.

!!! warning "Reflecting surfaces"
    It is important to consider that **all** surfaces of this mirror type are reflecting!
"""
struct Mirror{T, S <: AbstractShape{T}, C <: Tuple} <: AbstractReflectiveOptic{T}
    shape::S
    coatings::C
    function Mirror(shape::S, coatings::C = ()) where {T <: Real, S <: AbstractShape{T}, C <: Tuple}
        return new{T, S, C}(shape, coatings)
    end
end

"""
    SquarePlanoMirror2D(edge_length)

Constructs a 2D square plano [`Mirror`](@ref) with a given `edge_length`.
The reflecting surface is normal to the y-axis.

# Inputs

- `edge_length`: the edge length of the square mirror in [m]
"""
function SquarePlanoMirror2D(size::Real)
    T = float(typeof(size))
    shape = QuadraticFlatMesh(T(size))
    return Mirror(shape)
end

"""
    RectangularPlanoMirror(width, height, thickness)

Constructs a rectangular plano [`Mirror`](@ref) based on the input dimensions.
The front reflecting surface is normal to the y-axis and lies at the origin.

# Inputs

- `width`:      of the mirror in x-direction [m]
- `height`:     of the mirror in z-direction [m]
- `thickness`:  of the mirror in y-direction [m]
"""
function RectangularPlanoMirror(width::Real, height::Real, thickness::Real)
    T = float(promote_type(typeof(width), typeof(height), typeof(thickness)))
    shape = CuboidMesh(T(width), T(thickness), T(height))
    translate3d!(shape, [
        -T(width) / 2,    # x
        T(0),             # y
        -T(height) / 2    # z
    ])
    set_new_origin3d!(shape)
    return Mirror(shape)
end

"""
    SquarePlanoMirror(width, thickness)

Constructs a square plano [`Mirror`](@ref) with equal width and height.
The front reflecting surface is normal to the y-axis and lies at the origin.
See also [`RectangularPlanoMirror`](@ref).

# Inputs

- `width`: the side length of the square mirror in x- and y-direction [m]
- `thickness`: of the mirror in [m]
"""
function SquarePlanoMirror(width::Real, thickness::Real)
    return RectangularPlanoMirror(width, width, thickness)
end

"""
    RoundPlanoMirror <: AbstractReflectiveOptic

An ideal cylindrical mirror with planar reflecting surface, e.g. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`PlanoSurfaceSDF`](@ref) that represents the substrate
"""
struct RoundPlanoMirror{T} <: AbstractReflectiveOptic{T}
    shape::PlanoSurfaceSDF{T}
end

"""
    RoundPlanoMirror(diameter, thickness)

Returns a cylindrical, flat [`RoundPlanoMirror`](@ref) with perfect reflectivity based on:

# Inputs

- `diameter`: mirror diameter in [m]
- `thickness`: mirror substrate thickness in [m]
"""
function RoundPlanoMirror(diameter::Real, thickness::Real)
    T = float(promote_type(typeof(diameter), typeof(thickness)))
    shape = PlanoSurfaceSDF(T(thickness), T(diameter))
    return RoundPlanoMirror(shape)
end

"""[`ConcaveSphericalMirror`](@ref) shape type based on a [`UnionSDF`](@ref)"""
const ConcaveSphericalMirrorShape{T} = UnionSDF{
    T, Tuple{ConcaveSphericalSurfaceSDF{T}, PlanoSurfaceSDF{T}}}

"""
    ConcaveSphericalMirror <: AbstractReflectiveOptic

An ideal concave mirror with spherical reflecting surface, e.g. R = 1.
See also [`RoundPlanoMirror`](@ref).

# Fields

- `shape`: a [`ConcaveSphericalMirrorShape`](@ref) that represents the substrate
"""
struct ConcaveSphericalMirror{T} <: AbstractReflectiveOptic{T}
    shape::ConcaveSphericalMirrorShape{T}
end

"""
    ConcaveSphericalMirror(radius, thickness, diameter)

Constructor for a spherical mirror with a concave reflecting surface. The component is aligned with the positive y-axis.
See also [`ConcaveSphericalMirror`](@ref).

# Inputs

- `radius`: the spherical surface radius of curvature in [m]
- `thickness`: substrate thickness in [m]
- `diameter`: mirror outer diameter in [m]
"""
function ConcaveSphericalMirror(radius::Real, thickness::Real, diameter::Real)
    T = float(promote_type(typeof(radius), typeof(thickness), typeof(diameter)))
    cylinder = PlanoSurfaceSDF(T(thickness), T(diameter))
    concave = ConcaveSphericalSurfaceSDF(abs(T(radius)), T(diameter))
    shape = concave + cylinder
    return ConcaveSphericalMirror(shape)
end

"""
    RightAnglePrismMirror <: AbstractReflectiveOptic

An ideal right angle prism mirror with planar reflecting surface, i.e. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`RightAnglePrismSDF`](@ref) that represents the substrate
"""
struct RightAnglePrismMirror{T} <: AbstractReflectiveOptic{T}
    shape::RightAnglePrismSDF{T}
end

"""
    RightAnglePrismMirror(leg_length, height)

Constructs a right angle prism mirror. The primary surface is aligned with the pos. y-axis.

# Inputs

- `leg_length`: edge length in x and y in [m]
- `height`: in z-axis in [m]
"""
function RightAnglePrismMirror(leg_length::Real, height::Real)
    T = float(promote_type(typeof(leg_length), typeof(height)))
    shape = RightAnglePrismSDF(T(leg_length), T(height))
    zrotate3d!(shape, T(deg2rad(45 + 180)))
    return RightAnglePrismMirror(shape)
end
