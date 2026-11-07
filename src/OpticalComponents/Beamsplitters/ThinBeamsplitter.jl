"""
    ThinBeamsplitter <: AbstractBeamsplitter

Represents a 2D beam-splitting device: a thin layer whose [`BeamsplitterCoating`](@ref) acts at the
boundary on which it lies, see [`surface_model`](@ref).

# Fields

- `shape`: 2D [`AbstractShape`](@ref) at which the splitting process occurs (e.g. a 2D-[`Mesh`](@ref))
- `coating`: the [`BeamsplitterCoating`](@ref) with the field factors `reflectance` and `transmittance`

!!! warning
    The stored `reflectance` and `transmittance` are amplitudes. They should be calculated from the
    power reflectance ``R`` such that `reflectance^2 + transmittance^2 = 1`, as checked by `isvalid`.
"""
struct ThinBeamsplitter{T <: Real, S <: AbstractShape{T}} <: AbstractBeamsplitter{T}
    shape::S
    coating::BeamsplitterCoating{T}
end

# The factors in the number type of the `shape`, e.g. of a shape built from dual numbers for
# automatic differentiation w.r.t. its pose, while the reflectance is a `Float64`
function ThinBeamsplitter(shape::AbstractShape{T}, reflectance::Real, transmittance::Real) where {T}
    return ThinBeamsplitter{T, typeof(shape)}(shape, BeamsplitterCoating(T(reflectance), T(transmittance)))
end

surface_model(bs::ThinBeamsplitter) = bs.coating
reflectance(bs::ThinBeamsplitter) = reflectance(surface_model(bs))
transmittance(bs::ThinBeamsplitter) = transmittance(surface_model(bs))

"""
    ThinBeamsplitter(width, height; reflectance=0.5)

Creates a zero-thickness, lossless, non-polarizing 2D rectangular [`ThinBeamsplitter`](@ref) where

- `width`: is the x-dir. edge length in [m]
- `height`: is the z-dir. edge length in [m]
- `reflectance`: kw-arg in the open interval (0, 1) that determines how much light is **reflected**, i.e. 0.7 for a 70:30 splitter

# Additional information

!!! info "Reflectance"
    The input `reflectance` is the power reflectance ``R``. The amplitudes ``\\sqrt{R}`` and ``\\sqrt{1 - R}``
    are stored in the fields `reflectance` and `transmittance`, such that their squares sum to one.

!!! warning "Reflection phase jump"
    The reflection phase jump θᵣ = π applies to a beam that faces the normal of the splitter, see
    [`channels`](@ref) of its [`BeamsplitterCoating`](@ref).
"""
function ThinBeamsplitter(width::Real, height::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = RectangularFlatMesh(width, height)
    Reflected = sqrt(reflectance)
    Transmitted = sqrt(1 - Reflected^2)
    return ThinBeamsplitter(shape, Reflected, Transmitted)
end

"""
    RoundThinBeamsplitter(diameter; reflectance=0.5)

Creates a zero-thickness, 2D round [`ThinBeamsplitter`](@ref) with the specified `diameter` in [m].
For more information, refer to the [`ThinBeamsplitter`](@ref) constructor.
"""
function RoundThinBeamsplitter(diameter::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = CircularFlatMesh(diameter / 2)
    R = sqrt(reflectance)
    T = sqrt(1 - R^2)
    return ThinBeamsplitter(shape, R, T)
end

function ThinBeamsplitter(width::Real; reflectance::Real = 0.5)
    ThinBeamsplitter(width, width; reflectance)
end

Base.isvalid(bs::AbstractBeamsplitter) = reflectance(bs)^2 + transmittance(bs)^2 ≈ 1
