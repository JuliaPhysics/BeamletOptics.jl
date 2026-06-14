abstract type AbstractPlateComponent{T} <: AbstractObjectGroup{T} end

coating(p::AbstractPlateComponent) = p.coating
substrate(p::AbstractPlateComponent) = p.substrate

Base.position(p::AbstractPlateComponent) = position(coating(p))
orientation(p::AbstractPlateComponent) = orientation(substrate(p))
shape_trait_of(::AbstractPlateComponent) = MultiShape()
shape(p::AbstractPlateComponent) = (substrate(p), coating(p))
objects(p::AbstractPlateComponent) = (substrate(p), coating(p))
refractive_index(p::AbstractPlateComponent, λ::Real) = refractive_index(substrate(p), λ)

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
abstract type AbstractPlateBeamsplitter{T} <: AbstractPlateComponent{T} end

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
struct RectangularPlateBeamsplitter{T, M} <: AbstractPlateBeamsplitter{T}
    substrate::Prism{T, BoxSDF{T}}
    coating::Coating{T, Mesh{T}, M}
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
    T = float(promote_type(typeof(width), typeof(height), typeof(thickness)))
    # create substrate prism and move into pos
    substrate_shape = BoxSDF(T(width), T(thickness), T(height))
    substrate = Prism(substrate_shape, n)
    translate3d!(substrate, [T(0), T(thickness/2), T(0)])
    # rotate splitter "coating" into pos
    coating = ThinBeamsplitter(T(width), T(height); reflectance=T(reflectance))
    zrotate3d!(coating, T(π))
    M = typeof(coating.model)
    return RectangularPlateBeamsplitter{T, M}(substrate, coating)
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
struct RoundPlateBeamsplitter{T, M} <: AbstractPlateBeamsplitter{T}
    substrate::Prism{T, PlanoSurfaceSDF{T}}
    coating::Coating{T, Mesh{T}, M}
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
    T = float(promote_type(typeof(diameter), typeof(thickness)))
    # create substrate cylinder prism
    substrate_shape = PlanoSurfaceSDF(T(thickness), T(diameter))
    substrate = Prism(substrate_shape, n)
    # round splitter coating
    coating = RoundThinBeamsplitter(T(diameter); reflectance=T(reflectance))
    M = typeof(coating.model)
    return RoundPlateBeamsplitter{T, M}(substrate, coating)
end
