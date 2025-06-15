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
abstract type AbstractPlateBeamsplitter{T, S} <: AbstractBeamsplitter{T, S} end

coating(pbs::AbstractPlateBeamsplitter) = pbs.coating
substrate(pbs::AbstractPlateBeamsplitter) = pbs.substrate

Base.position(pbs::AbstractPlateBeamsplitter) = position(coating(pbs))

orientation(pbs::AbstractPlateBeamsplitter) = orientation(substrate(pbs))

shape_trait_of(::AbstractPlateBeamsplitter) = MultiShape()

shape(pbs::AbstractPlateBeamsplitter) = (substrate(pbs), coating(pbs))

refractive_index(pbs::AbstractPlateBeamsplitter, λ::Real) = refractive_index(substrate(pbs), λ)

"""Placeholder type for the  shape of a [`RectangularPlateBeamsplitter`](@ref)"""
struct RectangularPlateBeamsplitterShape{T} <: AbstractShape{T} end

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
struct RectangularPlateBeamsplitter{T} <: AbstractPlateBeamsplitter{T, RectangularPlateBeamsplitterShape{T}}
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

"""Placeholder type for the  shape of a [`RoundPlateBeamsplitterShape`](@ref)"""
struct RoundPlateBeamsplitterShape{T} <: AbstractShape{T} end

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
struct RoundPlateBeamsplitter{T} <: AbstractPlateBeamsplitter{T, RoundPlateBeamsplitterShape{T}}
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
    # this is sooooooo stupid but necessary to ensure correct intersection... 
    ic = intersect3d(coating(pbs), ray)
    is = intersect3d(substrate(pbs), ray)
    if isnothing(ic) & isnothing(is)
        return nothing
    end
    if isnothing(is)
        object!(ic, pbs)
        return ic
    end
    if isnothing(ic)
        object!(is, pbs)
        return is
    end
    # if both shapes are hit, pick the coating
    if length(ic) ≈ length(is)
        object!(ic, pbs)
        return ic
    end
    if length(ic) < length(is)
        object!(ic, pbs)
        return ic
    else
        object!(is, pbs)
        return is
    end
end

function interact3d(
    system::AbstractSystem,
    pbs::AbstractPlateBeamsplitter,
    beam::Beam{T, R},
    ray::R
    ) where {T <: Real, R <: AbstractRay{T}}
    # Substrate interaction
    if shape(intersection(ray)) === shape(substrate(pbs))
        interaction = interact3d(system, substrate(pbs), beam, ray)
        # Hint towards coating
        hint!(interaction, Hint(pbs, shape(coating(pbs))))
        return interaction
    end
    # Splitter "coating" interaction
    if shape(intersection(ray)) === shape(coating(pbs))
        # Beamsplitter coating interaction
        interact3d(system, coating(pbs), beam, ray)
        # Update refractive index and calculate refraction
        λ = wavelength(ray)
        n_optics = refractive_index(pbs, λ)
        n_system = refractive_index(system, λ)
        if isentering(ray)
            # transmitted ray is refracted into substrate
            _nt = n_optics
            _nr = n_system
            n_d, _ = refraction3d(ray, n_optics)
        else
            # transmitted ray is refracted into environment
            _nt = n_system
            _nr = n_optics
            n_d, _ = refraction3d(ray, n_system)
        end
        refractive_index!(first(rays(beam.children[1])), _nt)
        refractive_index!(first(rays(beam.children[2])), _nr)
        direction!(first(rays(beam.children[1])), n_d)
        return nothing
    end
    # if nothing worked, return nothing
    return nothing
end

function interact3d(
    system::AbstractSystem,
    pbs::AbstractPlateBeamsplitter,
    gauss::GaussianBeamlet,
    id::Int
    )
    _shape = shape(intersection(rays(gauss.chief)[id]))
    # Substrate interaction
    if _shape === shape(substrate(pbs))
        interaction = interact3d(system, substrate(pbs), gauss, id)
        hint!(interaction, Hint(pbs, shape(coating(pbs))))
        return interaction
    end
    # Splitter interaction
    if _shape === shape(coating(pbs))
        interact3d(system, coating(pbs), gauss, id)
        # Update refractive index and calculate refraction
        λ = wavelength(rays(gauss.chief)[id])
        n_optics = refractive_index(pbs, λ)
        n_system = refractive_index(system, λ)
        if isentering(gauss, id)
            # transmitted ray is refracted into substrate
            _nt = n_optics
            _nr = n_system
            n_c, _ = refraction3d(rays(gauss.chief)[id], n_optics) 
            n_w, _ = refraction3d(rays(gauss.waist)[id], n_optics) 
            n_d, _ = refraction3d(rays(gauss.divergence)[id], n_optics) 
        else
            # transmitted ray is refracted into environment
            _nt = n_system
            _nr = n_optics
            n_c, _ = refraction3d(rays(gauss.chief)[id], n_system) 
            n_w, _ = refraction3d(rays(gauss.waist)[id], n_system) 
            n_d, _ = refraction3d(rays(gauss.divergence)[id], n_system) 
        end
        # Update children ref. index and dir. due to refraction
        refractive_index!(gauss.children[1], 1, _nt)
        refractive_index!(gauss.children[2], 1, _nr)
        direction!(first(rays(gauss.children[1].chief)), n_c)
        direction!(first(rays(gauss.children[1].waist)), n_w)
        direction!(first(rays(gauss.children[1].divergence)), n_d)
        return nothing
    end
    # if nothing worked, return nothing
    return nothing
end