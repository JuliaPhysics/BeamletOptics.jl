"""Placeholder type for the  shape of a [`CubeBeamsplitter`](@ref)"""
struct CubeBeamsplitterShape{T} <: AbstractShape{T} end

"""
    CubeBeamsplitter <: AbstractBeamsplitter

A cuboid beamsplitter where the splitting interaction occurs between two [`RightAnglePrism`](@ref)s.
For more information refer to the [`AbstractPlateBeamsplitter`](@ref) docs.

# Fields

- `front`: the forward facing substrate, represented by a [`RightAnglePrism`](@ref)
- `back`: the backward facing substrate, represented by a [`RightAnglePrism`](@ref)
- `coating`: a rectangular [`ThinBeamsplitter`](@ref) that represents the splitting interface

# Additional information

!!! info "Hints and interaction logic"
    In order to model gap-free beam propagation, the `interact3d` model relies heavily on the [`Hint`](@ref)-API.
    If the `front` or `back` substrate is hit, the `Hint` will ensure that the beam intersects the `coating`.
"""
struct CubeBeamsplitter{T, N, C} <: AbstractBeamsplitter{T}
    front::Prism{T, RightAnglePrismSDF{T}, N, C}
    back::Prism{T, RightAnglePrismSDF{T}, N, C}
    coating::ThinBeamsplitter{T, Mesh{T}}
end

shape_trait_of(::CubeBeamsplitter) = MultiShape()

shape(cbs::CubeBeamsplitter) = (cbs.front, cbs.back, cbs.coating)

refractive_index(cbs::CubeBeamsplitter, λ::Real) = refractive_index(cbs.front, λ)

"""
    CubeBeamsplitter(leg_length, n; reflectance=0.5)

Creates a [`CubeBeamsplitter`](@ref). The cuboid is centered at the origin. The splitter
coating is orientated at a 45° angle with respect to the y-axis.

# Inputs

- `leg_length`: the x-, y- and z-edge length in [m]
- `n`: the [`RefractiveIndex`](@ref) of the front and back prism

# Keywords

- `reflectance`: defines the splitting ratio in [-], i.e. R = 0 ... 1.0
"""
function CubeBeamsplitter(
        leg_length::Real,
        n::RefractiveIndex;
        reflectance::Real = 0.5
)
    # Use AR coating for prism faces to ensure 100% transmission to the splitter
    ar = SimpleCoating(0.0)
    front = RightAnglePrism(leg_length, leg_length, n, ar)
    back = RightAnglePrism(leg_length, leg_length, n, ar)
    # Rotate back prism to form the other half of the cube
    zrotate3d!(back, π)

    bs = ThinBeamsplitter(√2 * leg_length, leg_length; reflectance)
    zrotate3d!(bs, deg2rad(180 - 45))
    set_new_origin3d!(shape(bs))

    return CubeBeamsplitter{typeof(leg_length), typeof(n), typeof(ar)}(front, back, bs)
end

function interact3d(
        system::AbstractSystem,
        cbs::CubeBeamsplitter,
        beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    # Front prism interaction
    if shape(intersection(ray)) === shape(cbs.front)
        interaction = interact3d(system, cbs.front, beam, ray)
        # Hint towards coating
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
    # Splitter "coating" interaction
    if shape(intersection(ray)) === shape(cbs.coating)
        # Beamsplitter coating interaction
        interact3d(system, cbs.coating, beam, ray)
        # Update refractive index
        _n = refractive_index(cbs, wavelength(ray))
        refractive_index!(first(rays(beam.children[1])), _n)
        refractive_index!(first(rays(beam.children[2])), _n)
        return nothing
    end
    # Back prism interaction
    if shape(intersection(ray)) === shape(cbs.back)
        interaction = interact3d(system, cbs.back, beam, ray)
        # Hint towards coating
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
end

function interact3d(
        system::AbstractSystem,
        cbs::CubeBeamsplitter,
        gauss::GaussianBeamlet,
        id::Int)
    _shape = shape(intersection(rays(gauss.chief)[id]))
    # Front prism interaction
    if _shape === shape(cbs.front)
        interaction = interact3d(system, cbs.front, gauss, id)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
    # Splitter "coating" interaction
    if _shape === shape(cbs.coating)
        interact3d(system, cbs.coating, gauss, id)
        # Update refractive index
        _n = refractive_index(cbs, wavelength(gauss))
        refractive_index!(gauss.children[1], 1, _n)
        refractive_index!(gauss.children[2], 1, _n)
        return nothing
    end
    # Back prism interaction
    if _shape === shape(cbs.back)
        interaction = interact3d(system, cbs.back, gauss, id)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
end
