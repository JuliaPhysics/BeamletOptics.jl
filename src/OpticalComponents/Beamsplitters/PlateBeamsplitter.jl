"""
    AbstractPlateComponent <: AbstractObject

Generic type for plate-like optical components.
"""
abstract type AbstractPlateComponent{T} <: AbstractObject{T} end

"""
    AbstractPlateBeamsplitter <: AbstractRefractiveOptic

A generic type to represent an optical plate beamsplitter where beam splitting occurs at a coated surface face.
"""
abstract type AbstractPlateBeamsplitter{T, N} <: AbstractRefractiveOptic{T, N} end


"""
    RectangularPlateBeamsplitter{T, S, N, C} <: AbstractPlateBeamsplitter{T, N}

A plate beamsplitter with a rectangular substrate and surface coatings.
"""
struct RectangularPlateBeamsplitter{T, S <: BoxSDF{T}, N <: RefractiveIndex, C <: Tuple} <: AbstractPlateBeamsplitter{T, N}
    shape::S
    n::N
    coatings::C
    function RectangularPlateBeamsplitter(
            shape::S, n::N, coatings::C = ()) where {T <: Real, S <: BoxSDF{T}, N <: RefractiveIndex, C <: Tuple}
        test_refractive_index_function(n)
        return new{T, S, N, C}(shape, n, coatings)
    end
end

"""
    RoundPlateBeamsplitter{T, S, N, C} <: AbstractPlateBeamsplitter{T, N}

A plate beamsplitter with a cylindrical substrate and surface coatings.
"""
struct RoundPlateBeamsplitter{T, S <: PlanoSurfaceSDF{T}, N <: RefractiveIndex, C <: Tuple} <: AbstractPlateBeamsplitter{T, N}
    shape::S
    n::N
    coatings::C
    function RoundPlateBeamsplitter(
            shape::S, n::N, coatings::C = ()) where {T <: Real, S <: PlanoSurfaceSDF{T}, N <: RefractiveIndex, C <: Tuple}
        test_refractive_index_function(n)
        return new{T, S, N, C}(shape, n, coatings)
    end
end

# Compatibility accessors
substrate(p::AbstractPlateBeamsplitter) = p
coating(p::AbstractPlateBeamsplitter) = get_matching_coating(coatings(p), shape(p), [0.0, 0.0, 0.0], [0.0, -1.0, 0.0])

"""
    RectangularPlateBeamsplitter(width, height, thickness, n; reflectance=0.5, back_coating=nothing)

Creates a [`RectangularPlateBeamsplitter`](@ref). The splitter front face is centered at the origin,
with the substrate extending along the positive y-axis.

# Inputs

- `width`: substrate width along the x-axis in [m]
- `height`: substrate height along the z-axis in [m]
- `thickness`: substrate thickness along the y-axis in [m]
- `n`: the [`RefractiveIndex`](@ref) of the substrate

# Keywords

- `reflectance`: defines the splitting ratio in [-], i.e. R = 0 ... 1.0
- `back_coating`: optional coating attached to the back face (e.g., [`SimpleARCoating`](@ref))
"""
function RectangularPlateBeamsplitter(
        width::Real,
        height::Real,
        thickness::Real,
        n::RefractiveIndex;
        reflectance::Real=0.5,
        back_coating=nothing
    )
    if reflectance >= 1 || reflectance <= 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    T = float(promote_type(typeof(width), typeof(height), typeof(thickness)))
    substrate_shape = BoxSDF(T(width), T(thickness), T(height))
    translate3d!(substrate_shape, [T(0), T(thickness/2), T(0)])

    r = sqrt(reflectance)
    t = sqrt(1.0 - reflectance)
    coat_bs = SimpleBeamsplitterCoating(r, r, t, t)

    coatings_list = if isnothing(back_coating)
        (:front => coat_bs,)
    else
        (:front => coat_bs, :back => _coating_model(back_coating))
    end

    return RectangularPlateBeamsplitter(substrate_shape, n, coatings_list)
end

"""
    RoundPlateBeamsplitter(diameter, thickness, n; reflectance=0.5, back_coating=nothing)

Creates a circular [`RoundPlateBeamsplitter`](@ref). The splitter front face is centered at the origin,
with the substrate extending along the positive y-axis.

# Inputs

- `diameter`: x-z-plane substrate diameter in [m]
- `thickness`: substrate thickness along the y-axis in [m]
- `n`: the [`RefractiveIndex`](@ref) of the substrate

# Keywords

- `reflectance`: defines the splitting ratio in [-], i.e. R = 0 ... 1.0
- `back_coating`: optional coating attached to the back face (e.g., [`SimpleARCoating`](@ref))
"""
function RoundPlateBeamsplitter(
        diameter::Real,
        thickness::Real,
        n::RefractiveIndex;
        reflectance::Real=0.5,
        back_coating=nothing
    )
    if reflectance >= 1 || reflectance <= 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    T = float(promote_type(typeof(diameter), typeof(thickness)))
    substrate_shape = PlanoSurfaceSDF(T(thickness), T(diameter))

    r = sqrt(reflectance)
    t = sqrt(1.0 - reflectance)
    coat_bs = SimpleBeamsplitterCoating(r, r, t, t)

    coatings_list = if isnothing(back_coating)
        (:front => coat_bs,)
    else
        (:front => coat_bs, :back => _coating_model(back_coating))
    end

    return RoundPlateBeamsplitter(substrate_shape, n, coatings_list)
end


