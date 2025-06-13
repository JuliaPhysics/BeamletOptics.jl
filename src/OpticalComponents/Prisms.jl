"""
    Prism{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, S, N}

Essentially represents the same functionality as [`Lens`](@ref).
Refer to its documentation.
"""
struct Prism{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, S, N}
    shape::S
    n::N
    function Prism(shape::S, n::N) where {T<:Real, S<:AbstractShape{T}, N<:RefractiveIndex}
        test_refractive_index_function(n)
        return new{T, S, N}(shape, n)
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
function RightAnglePrism(leg_length::Real, height::Real, n::BeamletOptics.RefractiveIndex)
    shape = BeamletOptics.RightAnglePrismSDF(leg_length, height)
    return Prism(shape, n)
end
