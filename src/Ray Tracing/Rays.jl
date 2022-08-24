mutable struct Intersection{T}
    const t::T
    const n::Union{Vector{T}, Missing}
    oID::Union{Int, Missing}
end

# Optimized for pointer look-up
const _NoIntersectionF64 = Intersection(Inf64, missing, missing)
const _NoIntersectionF32 = Intersection(Inf32, missing, missing)

NoIntersection(::Type{Float64}) = _NoIntersectionF64
NoIntersection(::Type{Float32}) = _NoIntersectionF32

"""
    Ray{T}

Mutable struct to store ray information. A `Ray` is described by ``\\vec{v}_{pos}+t\\cdot\\vec{v}_{dir}`` with ``t\\in[0,\\infty)``.

# Fields
- `pos`: a 3D-vector that describes the `Ray` origin
- `dir`: a normalized 3D-vector that describes the `Ray` direction
- `intersection`: refer to `Intersection{T}`.
"""
mutable struct Ray{T}
    pos::Vector{T}
    dir::Vector{T}
    intersection::Intersection{T}
end

"""
    Ray(pos::Vector{P}, dir::Vector{D}) where {P, D}

Constructs an instance of `Ray` with common type `T`. Vector `dir` is normalized. Ray is initialized with `NoIntersection`.\\
In general, `T<:AbstractFloat`! Inputs of type `Int` are promoted to `Float64`.
"""
function Ray(pos::Vector{P}, dir::Vector{D}) where {P<:Real, D<:Real}
    @assert norm3d(dir) != 0 "Illegal vector for direction"
    T = promote_type(P, D)
    # Catch integer type inputs
    if !(T<:AbstractFloat)
        T = Float64
    end
    return Ray{T}(T.(pos), normalize3d(T.(dir)), NoIntersection(T))
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