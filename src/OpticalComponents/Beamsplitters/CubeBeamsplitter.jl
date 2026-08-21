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
struct CubeBeamsplitter{T} <: AbstractBeamsplitter{T}
    front::Prism{T, RightAnglePrismSDF{T}}
    back::Prism{T, RightAnglePrismSDF{T}}
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
        reflectance::Real=0.5
    )
    front = RightAnglePrism(leg_length, leg_length, n)
    back = RightAnglePrism(leg_length, leg_length, n)
    bs = ThinBeamsplitter(√2*leg_length, leg_length; reflectance)
    zrotate3d!(back, deg2rad(180))
    zrotate3d!(bs, deg2rad(180-45))
    set_new_origin3d!(shape(bs))
    return CubeBeamsplitter(front, back, bs)
end

function _cbs_hit_target(cbs::CubeBeamsplitter, ray::AbstractRay, int::Nullable{<:AbstractIntersection} = nothing)
    if isnothing(int)
        int = intersect3d(shape(cbs), ray)
    end
    isnothing(int) && return :none
    tol = Config.get_coincident_boundary_tolerance()
    ic = intersect3d(shape(cbs.coating), ray)
    if !isnothing(ic) && isapprox(length(int), length(ic), atol = tol)
        return :coating
    end
    if_hit = intersect3d(shape(cbs.front), ray)
    if !isnothing(if_hit) && isapprox(length(int), length(if_hit), atol = tol)
        return :front
    end
    ib = intersect3d(shape(cbs.back), ray)
    if !isnothing(ib) && isapprox(length(int), length(ib), atol = tol)
        return :back
    end
    return :none
end

function interact3d(
    system::AbstractSystem,
    cbs::CubeBeamsplitter,
    beam::Beam{T, R},
    ray::R) where {T <: Real, R <: AbstractRay{T}}
    int = isempty(beam.segments) ? intersect3d(shape(cbs), ray) : last(beam.segments).intersection
    target = _cbs_hit_target(cbs, ray, int)
    if target === :front
        interaction = interact3d(system, cbs.front, beam, ray)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    elseif target === :coating
        interact3d(system, cbs.coating, beam, ray)
        _n = refractive_index(cbs, wavelength(ray))
        for child in beam.children
            if child.head_ray isa Ray
                child.head_ray = Ray(position(child.head_ray), direction(child.head_ray), wavelength(child.head_ray), _n)
            elseif child.head_ray isa PolarizedRay
                child.head_ray = PolarizedRay(position(child.head_ray), direction(child.head_ray), wavelength(child.head_ray), _n, polarization(child.head_ray))
            end
        end
        return nothing
    elseif target === :back
        interaction = interact3d(system, cbs.back, beam, ray)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
    return nothing
end

function interact3d(
    system::AbstractSystem,
    cbs::CubeBeamsplitter,
    gauss::GaussianBeamlet,
    id::Int)
    chief_ray = rays(gauss.chief)[id]
    int = id <= length(gauss.chief.segments) ? gauss.chief.segments[id].intersection : intersect3d(shape(cbs), chief_ray)
    target = _cbs_hit_target(cbs, chief_ray, int)
    if target === :front
        interaction = interact3d(system, cbs.front, gauss, id)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    elseif target === :coating
        interact3d(system, cbs.coating, gauss, id)
        _n = refractive_index(cbs, wavelength(gauss))
        for child in gauss.children
            child.chief.head_ray = Ray(position(child.chief.head_ray), direction(child.chief.head_ray), wavelength(child.chief.head_ray), _n)
            child.waist.head_ray = Ray(position(child.waist.head_ray), direction(child.waist.head_ray), wavelength(child.waist.head_ray), _n)
            child.divergence.head_ray = Ray(position(child.divergence.head_ray), direction(child.divergence.head_ray), wavelength(child.divergence.head_ray), _n)
        end
        return nothing
    elseif target === :back
        interaction = interact3d(system, cbs.back, gauss, id)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
    return nothing
end

function interact3d(
    system::AbstractSystem,
    cbs::CubeBeamsplitter,
    agb::AstigmaticGaussianBeamlet,
    id::Int)
    chief_ray = rays(agb.c)[id]
    int = id <= length(agb.c.segments) ? agb.c.segments[id].intersection : intersect3d(shape(cbs), chief_ray)
    target = _cbs_hit_target(cbs, chief_ray, int)
    if target === :front
        interaction = interact3d(system, cbs.front, agb, id)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    elseif target === :coating
        interact3d(system, cbs.coating, agb, id)
        _n = refractive_index(cbs, wavelength(agb))
        for child in agb.children
            for b in _component_beams(child)
                if b.head_ray isa PolarizedRay
                    b.head_ray = PolarizedRay(position(b.head_ray), direction(b.head_ray), wavelength(b.head_ray), _n, polarization(b.head_ray))
                else
                    b.head_ray = Ray(position(b.head_ray), direction(b.head_ray), wavelength(b.head_ray), _n)
                end
            end
        end

        return nothing
    elseif target === :back
        interaction = interact3d(system, cbs.back, agb, id)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
    return nothing
end