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