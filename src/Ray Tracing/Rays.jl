"""
    Ray{T}

A mutable struct that contains a **position vector** `pos`, a **directional vector** `dir`
and a **length** variable `t` that are used to describe a generic ray as `pos+t*dir`.
The directional vector is required/adjusted to have unit length, i.e. `abs(dir) == 1`.
"""
mutable struct Ray{T}
    pos::Vector{T}
    dir::Vector{T}
    len::T
    @doc """
        Ray{T}(pos::Vector{T}, dir::Vector{T}) where T

    Parametric type constructor for struct Ray. Takes in a position vector `pos` and directional vector `dir`, which is scaled
    to unit length. The initial ray length `t` is set to `Inf`.
    """
    function Ray{T}(pos::Vector{T}, dir::Vector{T}) where {T}
        @assert norm(dir) != 0 "Illegal vector for direction"
        new{T}(pos, dir / norm(dir), Inf)
    end
end

"""
    Ray(pos::Vector{T}, dir::Vector{G}) where {T<:Union{Int, Float64}, G<:Union{Int, Float64}}

Concrete implementation of Ray struct for `Float64` data type. Accepts all combinations where `pos` and `dir`
are of type `Float64` and/or abstract `Int`. Allows seperate implementation for `Float32`, etc.
Promotes input types to `Float64`.
"""
function Ray(pos::Vector{T}, dir::Vector{G}) where {T<:Union{Int,Float64},G<:Union{Int,Float64}}
    return Ray{Float64}(Float64.(pos), Float64.(dir))
end

mutable struct Beamlet{T}
    chief::Vector{Ray{T}}
    divergence::Vector{Ray{T}}
    waist::Vector{Ray{T}}

    # Beamlet constructor
    # sizehint! in constructor?
    # append! performance (irrelevant)
    function Beamlet(chief::Ray, λ, w0; T=Float64)
        # All rays are initialized parallel to the x,y-plane
        # Divergence angle in rad
        θ = λ / (π * w0)
        # Divergence ray
        divergence = Ray(chief.pos, SCDI.rotate3d([0, 0, 1], θ) * chief.dir)
        # Waist ray
        waist = Ray(chief.pos + SCDI.orthogonal3d(chief.dir, [0, 0, 1]) * w0, chief.dir)
        new{T}([chief], [divergence], [waist])
    end
end

"""
    Beam

Temporary container struct to test ray tracing. Wavelength `λ` in nm.
"""
mutable struct Beam
    rays::Vector{Ray{Float64}}
    λ::Float64
end