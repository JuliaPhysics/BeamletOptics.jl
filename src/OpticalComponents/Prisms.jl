"""
    Prism{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, N}

Essentially represents the same functionality as [`Lens`](@ref).
Refer to its documentation.
"""
struct Prism{T, S <: AbstractShape{T}, N <: RefractiveIndex, C <: Tuple} <: AbstractRefractiveOptic{T, N}
    shape::S
    n::N
    coatings::C
    function Prism(shape::S, n::N, coatings::C = ()) where {T <: Real, S <: AbstractShape{T}, N <: RefractiveIndex, C <: Tuple}
        test_refractive_index_function(n)
        return new{T, S, N, C}(shape, n, coatings)
    end
end

thickness(p::Prism) = thickness(shape(p))

"""
    RightAnglePrism(leg_length, height, n)

Creates a right angle symmetric [`Prism`](@ref). The prism is *not aligned* with the y-axis.

# Inputs

- `leg_length`: dimension in x- and y-direction in [m]
- `height`: in [m]
- `n`: [`RefractiveIndex`](@ref) of the prism
"""
function RightAnglePrism(leg_length::Real, height::Real, n::RefractiveIndex)
    shape = RightAnglePrismSDF(leg_length, height)
    return Prism(shape, n)
end
