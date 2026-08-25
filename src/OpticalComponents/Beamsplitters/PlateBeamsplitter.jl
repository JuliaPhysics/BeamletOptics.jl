"""
    AbstractPlateBeamsplitter{T} <: AbstractObject{T}

A generic type to represent an [`AbstractBeamsplitter`](@ref) composite component that consists of a substrate with a
single coated face at which a beam splitting interaction occurs.

# Implementation reqs.

Subtypes of `AbstractPlateBeamsplitter` should implement all supertype reqs. as well as:

## Fields

- `coating`: a [`ThinBeamsplitter`](@ref) that represents the splitter coating
- `substrate`: a [`Prism`](@ref) that represents the substrate

## Getters/setters

- `coating`: returns a [`ThinBeamsplitter`](@ref)
- `substrate`: returns a [`Prism`](@ref)
"""
abstract type AbstractPlateBeamsplitter{T} <: AbstractObject{T} end

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
