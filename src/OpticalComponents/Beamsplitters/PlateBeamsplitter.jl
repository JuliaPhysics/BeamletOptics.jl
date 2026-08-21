"""
    AbstractPlateBeamsplitter <: AbstractBeamsplitter

A generic type to represent an [`AbstractBeamsplitter`](@ref) that consists of a substrate with a 
single coated face at which a beam splitting interaction occurs.

# Implementation reqs.

Subtypes of `AbstractPlateBeamsplitter` should implement all supertype reqs. as well as:

## Fields

- `coating`: a [`ThinBeamsplitter`](@ref) that represents the splitter coating
- `substrate`: a [`Prism`](@ref) that represents the substrate

## Getters/setters

If the concrete implementation does not define the above fields, the following getters must be defined:

- `coating`: returns a [`ThinBeamsplitter`](@ref)
- `substrate`: returns a [`Prism`](@ref)

# Additional information

!!! info "Object orientation"
    This `interact3d` method of this type strongly assumes that the coating is positioned directly upon
    a single face of the substrate with a 100% fill factor.

!!! info "Interaction logic"
    This type uses the [`Hint`](@ref)-API in order to ensure that the splitting interaction is correctly
    triggered at the coating.
"""
abstract type AbstractPlateBeamsplitter{T} <: AbstractBeamsplitter{T} end

coating(pbs::AbstractPlateBeamsplitter) = pbs.coating
substrate(pbs::AbstractPlateBeamsplitter) = pbs.substrate

Base.position(pbs::AbstractPlateBeamsplitter) = position(coating(pbs))

orientation(pbs::AbstractPlateBeamsplitter) = orientation(substrate(pbs))

shape_trait_of(::AbstractPlateBeamsplitter) = MultiShape()

shape(pbs::AbstractPlateBeamsplitter) = (substrate(pbs), coating(pbs))

refractive_index(pbs::AbstractPlateBeamsplitter, λ::Real) = refractive_index(substrate(pbs), λ)

"""
    RectangularPlateBeamsplitter <: AbstractPlateBeamsplitter

A plate beamsplitter with rectangular substrate and a single coated face.
For more information refer to the [`AbstractPlateBeamsplitter`](@ref) docs.

# Fields

- `substrate`: a rectangular [`Prism`](@ref) that acts as the substrate
- `coating`: a [`ThinBeamsplitter`](@ref) that acts as the coating

# Additional information

!!! info "Kinematic center"
    The center of kinematics of this splitter lies at the center of the coating.
"""
struct RectangularPlateBeamsplitter{T} <: AbstractPlateBeamsplitter{T}
    substrate::Prism{T, BoxSDF{T}}
    coating::ThinBeamsplitter{T, Mesh{T}}
end

"""
    RectangularPlateBeamsplitter(width, height, thickness, n; reflectance=0.5)

Creates a [`RectangularPlateBeamsplitter`](@ref). The splitter is aligned with the negative y-axis.
The splitter coating is centered at the origin. See also [`RoundPlateBeamsplitter`](@ref).

# Inputs

- `width`: substrate width along the x-axis in [m]
- `height`: substrate height along the z-axis in [m]
- `thickness`: substrate thickness along the y-axis in [m]
- `n`: the [`RefractiveIndex`](@ref) of the substrate

# Keywords

- `reflectance`: defines the splitting ratio in [-], i.e. R = 0 ... 1.0
"""
function RectangularPlateBeamsplitter(
        width::Real,
        height::Real,
        thickness::Real,
        n::RefractiveIndex;
        reflectance::Real=0.5
    )
    # create substrate prism and move into pos
    substrate_shape = BoxSDF(width, thickness, height)
    substrate = Prism(substrate_shape, n)
    translate3d!(substrate, [0, thickness/2, 0])
    # rotate splitter "coating" into pos
    coating = ThinBeamsplitter(width, height; reflectance)
    zrotate3d!(coating, π)
    return RectangularPlateBeamsplitter(substrate, coating)
end

"""
    RoundPlateBeamsplitter <: AbstractPlateBeamsplitter

A plate beamsplitter with cylindrical substrate and a single coated face.
For more information refer to the [`AbstractPlateBeamsplitter`](@ref) docs.

# Fields

- `substrate`: a cylindrical [`Prism`](@ref) that acts as the substrate
- `coating`: a [`RoundThinBeamsplitter`](@ref) that acts as the coating

# Additional information

!!! info "Kinematic center"
    The center of kinematics of this splitter lies at the center of the coating.
"""
struct RoundPlateBeamsplitter{T} <: AbstractPlateBeamsplitter{T}
    substrate::Prism{T, PlanoSurfaceSDF{T}}
    coating::ThinBeamsplitter{T, Mesh{T}}
end

"""
    RoundPlateBeamsplitter(diameter, thickness, n; reflectance=0.5)

Creates a [`RoundPlateBeamsplitter`](@ref). The splitter is aligned with the negative y-axis.
The coating is centered at the origin. See also [`RectangularPlateBeamsplitter`](@ref).

# Inputs

- `diameter`: x-z-plane substrate diameter in [m]
- `thickness`: substrate thickness along the z-axis in [m]
- `n`: the [`RefractiveIndex`](@ref) of the substrate

# Keywords 

- `reflectance`: defines the splitting ratio in [-], i.e. R = 0 ... 1.0
"""
function RoundPlateBeamsplitter(
        diameter::Real,
        thickness::Real,
        n::RefractiveIndex;
        reflectance::Real=0.5
    )
    # create substrate cylinder prism
    substrate_shape = PlanoSurfaceSDF(thickness, diameter)
    substrate = Prism(substrate_shape, n)
    # round splitter coating (neg. y-axi normals)
    coating = RoundThinBeamsplitter(diameter; reflectance)
    return RoundPlateBeamsplitter(substrate, coating)
end

function intersect3d(pbs::AbstractPlateBeamsplitter, ray::AbstractRay)
    ic = intersect3d(coating(pbs), ray)
    is = intersect3d(substrate(pbs), ray)
    if isnothing(ic) && isnothing(is)
        return nothing
    elseif isnothing(is)
        return ic
    elseif isnothing(ic)
        return is
    elseif length(ic) <= length(is)
        return ic
    else
        return is
    end
end

function _is_coating_hit(pbs::AbstractPlateBeamsplitter, ray::AbstractRay, int::Nullable{<:AbstractIntersection} = nothing)
    if isnothing(int)
        int = intersect3d(shape(pbs), ray)
    end
    isnothing(int) && return false
    ic = intersect3d(coating(pbs), ray)
    tol = Config.get_coincident_boundary_tolerance()
    return !isnothing(ic) && isapprox(length(int), length(ic), atol = tol)
end

function interact3d(
    system::AbstractSystem,
    pbs::AbstractPlateBeamsplitter,
    beam::Beam{T, R},
    ray::R
    ) where {T <: Real, R <: AbstractRay{T}}
    int = isempty(beam.segments) ? intersect3d(shape(pbs), ray) : last(beam.segments).intersection
    hit_coating = _is_coating_hit(pbs, ray, int)

    # Substrate interaction
    if !hit_coating
        interaction = interact3d(system, substrate(pbs), beam, ray)
        hint!(interaction, Hint(pbs, shape(coating(pbs))))
        return interaction
    end
    # Splitter "coating" interaction
    if hit_coating
        # Beamsplitter coating interaction
        interact3d(system, coating(pbs), beam, ray)
        # Update refractive index and calculate refraction
        λ = wavelength(ray)
        n_optics = refractive_index(pbs, λ)
        n_system = refractive_index(system, λ)
        nml = normal3d(int)
        entering = isentering(direction(ray), nml)
        if entering
            # transmitted ray is refracted into substrate
            _nt = n_optics
            _nr = n_system
            n_d, _ = refraction3d(direction(ray), nml, refractive_index(ray), n_optics)
        else
            # transmitted ray is refracted into environment
            _nt = n_system
            _nr = n_optics
            n_d, _ = refraction3d(direction(ray), -nml, refractive_index(ray), n_system)
        end
        child_t = beam.children[1]
        child_r = beam.children[2]
        if child_t.head_ray isa Ray
            child_t.head_ray = Ray(position(child_t.head_ray), n_d, wavelength(child_t.head_ray), _nt)
            child_r.head_ray = Ray(position(child_r.head_ray), direction(child_r.head_ray), wavelength(child_r.head_ray), _nr)
        elseif child_t.head_ray isa PolarizedRay
            child_t.head_ray = PolarizedRay(position(child_t.head_ray), n_d, wavelength(child_t.head_ray), _nt, polarization(child_t.head_ray))
            child_r.head_ray = PolarizedRay(position(child_r.head_ray), direction(child_r.head_ray), wavelength(child_r.head_ray), _nr, polarization(child_r.head_ray))
        end
        return nothing
    end
    return nothing
end

function interact3d(
    system::AbstractSystem,
    pbs::AbstractPlateBeamsplitter,
    gauss::GaussianBeamlet,
    id::Int
    )
    chief_ray = rays(gauss.chief)[id]
    int = id <= length(gauss.chief.segments) ? gauss.chief.segments[id].intersection : intersect3d(shape(pbs), chief_ray)
    hit_coating = _is_coating_hit(pbs, chief_ray, int)

    # Substrate interaction
    if !hit_coating
        interaction = interact3d(system, substrate(pbs), gauss, id)
        hint!(interaction, Hint(pbs, shape(coating(pbs))))
        return interaction
    end
    # Splitter interaction
    if hit_coating
        interact3d(system, coating(pbs), gauss, id)
        # Update refractive index and calculate refraction
        λ = wavelength(chief_ray)
        n_optics = refractive_index(pbs, λ)
        n_system = refractive_index(system, λ)
        nml = normal3d(int)
        entering = isentering(direction(chief_ray), nml)
        if entering
            # transmitted ray is refracted into substrate
            _nt = n_optics
            _nr = n_system
            n_c, _ = refraction3d(direction(rays(gauss.chief)[id]), nml, refractive_index(rays(gauss.chief)[id]), n_optics) 
            n_w, _ = refraction3d(direction(rays(gauss.waist)[id]), nml, refractive_index(rays(gauss.waist)[id]), n_optics) 
            n_d, _ = refraction3d(direction(rays(gauss.divergence)[id]), nml, refractive_index(rays(gauss.divergence)[id]), n_optics) 
        else
            # transmitted ray is refracted into environment
            _nt = n_system
            _nr = n_optics
            n_c, _ = refraction3d(direction(rays(gauss.chief)[id]), -nml, refractive_index(rays(gauss.chief)[id]), n_system) 
            n_w, _ = refraction3d(direction(rays(gauss.waist)[id]), -nml, refractive_index(rays(gauss.waist)[id]), n_system) 
            n_d, _ = refraction3d(direction(rays(gauss.divergence)[id]), -nml, refractive_index(rays(gauss.divergence)[id]), n_system) 
        end
        # Update children head rays
        c1 = gauss.children[1]
        c2 = gauss.children[2]
        c1.chief.head_ray = Ray(position(c1.chief.head_ray), n_c, wavelength(c1.chief.head_ray), _nt)
        c1.waist.head_ray = Ray(position(c1.waist.head_ray), n_w, wavelength(c1.waist.head_ray), _nt)
        c1.divergence.head_ray = Ray(position(c1.divergence.head_ray), n_d, wavelength(c1.divergence.head_ray), _nt)
        c2.chief.head_ray = Ray(position(c2.chief.head_ray), direction(c2.chief.head_ray), wavelength(c2.chief.head_ray), _nr)
        c2.waist.head_ray = Ray(position(c2.waist.head_ray), direction(c2.waist.head_ray), wavelength(c2.waist.head_ray), _nr)
        c2.divergence.head_ray = Ray(position(c2.divergence.head_ray), direction(c2.divergence.head_ray), wavelength(c2.divergence.head_ray), _nr)
        return nothing
    end
    return nothing
end

function interact3d(
    system::AbstractSystem,
    pbs::AbstractPlateBeamsplitter,
    agb::AstigmaticGaussianBeamlet,
    id::Int
    )
    chief_ray = rays(agb.c)[id]
    int = id <= length(agb.c.segments) ? agb.c.segments[id].intersection : intersect3d(shape(pbs), chief_ray)
    hit_coating = _is_coating_hit(pbs, chief_ray, int)

    # Substrate interaction
    if !hit_coating
        interaction = interact3d(system, substrate(pbs), agb, id)
        hint!(interaction, Hint(pbs, shape(coating(pbs))))
        return interaction
    end
    # Splitter interaction
    if hit_coating
        interact3d(system, coating(pbs), agb, id)
        # Update refractive index and calculate refraction
        λ = wavelength(chief_ray)
        n_optics = refractive_index(pbs, λ)
        n_system = refractive_index(system, λ)
        nml = normal3d(int)
        entering = isentering(direction(chief_ray), nml)
        if entering
            _nt = n_optics
            _nr = n_system
            n_target = n_optics
            n_normal = nml
        else
            _nt = n_system
            _nr = n_optics
            n_target = n_system
            n_normal = -nml
        end
        # Update children head rays
        for (beam, pbeam) in zip(_component_beams(agb.children[1]), _component_beams(agb))
            p_ray = rays(pbeam)[id]
            n_d, _ = refraction3d(direction(p_ray), n_normal, refractive_index(p_ray), n_target)
            if beam.head_ray isa PolarizedRay
                beam.head_ray = PolarizedRay(position(beam.head_ray), n_d, wavelength(beam.head_ray), _nt, polarization(beam.head_ray))
            else
                beam.head_ray = Ray(position(beam.head_ray), n_d, wavelength(beam.head_ray), _nt)
            end
        end
        for beam in _component_beams(agb.children[2])
            if beam.head_ray isa PolarizedRay
                beam.head_ray = PolarizedRay(position(beam.head_ray), direction(beam.head_ray), wavelength(beam.head_ray), _nr, polarization(beam.head_ray))
            else
                beam.head_ray = Ray(position(beam.head_ray), direction(beam.head_ray), wavelength(beam.head_ray), _nr)
            end
        end
        return nothing

    end
    return nothing
end

