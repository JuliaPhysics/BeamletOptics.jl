"""
    Intersection{T}

Stores some data resulting from ray tracing a `System`. This information can be used, i.e. for retracing.

# Fields:
- `t`: length of the ray parametrization in [m]
- `n`: normal vector at the point of intersection
- `oID`: index of the intersected object in the `System` object-vector
"""
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
    Information{T}

Stores the optical parameters that are relevant for the propagation of a ray through an optical system. See also `Ray{T}`.

# Fields:
- `λ`: wavelength in [nm]
- `n`: refractive index along the beam path
- `I`: intensity value [**currently a placeholder**]
- `P`: polarization vector (either Stokes or Jones formalism, tbd.) [**currently a placeholder**]
"""
mutable struct Information{T}
    λ::Union{T, Missing}
    n::Union{T, Missing}
    I::Union{T, Missing}
    P::Union{Vector{Complex{T}}, Missing}
end

# Optimized for pointer look-up
const _NoInformationF64 = Information{Float64}(missing, missing, missing, missing)
const _NoInformationF32 = Information{Float32}(missing, missing, missing, missing)

NoInformation(::Type{Float64}) = _NoInformationF64
NoInformation(::Type{Float32}) = _NoInformationF32

"Prototype constructor for `Information{T}`"
function Information(λ::T, n=1, I=1, P=[0,0]) where T
    return Information{T}(λ, n, I, P)
end

"""
    Ray{T}

Mutable struct to store ray information. A `Ray` is described by ``\\vec{v}_{pos}+t\\cdot\\vec{v}_{dir}`` with ``t\\in[0,\\infty)``.

# Fields
- `pos`: a 3D-vector that describes the `Ray` origin
- `dir`: a normalized 3D-vector that describes the `Ray` direction
- `intersection`: refer to `Intersection{T}`.
- `information`: refer to `Information{T}`
"""
mutable struct Ray{T} <: AbstractRay{T}
    pos::Vector{T}
    dir::Vector{T}
    intersection::Intersection{T}
    information::Information{T}
end

Base.length(ray::Ray) = ray.intersection.t
length!(ray::Ray, len) = nothing

information(ray::Ray) = ray.information
information!(ray::Ray, info) = (ray.information = info)

wavelength(ray::Ray) = ray.information.λ
wavelength!(ray::Ray, λ) = (ray.information.λ = λ)

refractive_index(ray::Ray) = ray.information.n
refractive_index!(ray::Ray, n) = (ray.information.n = n)

# Placeholders for the future
polarization(ray::Ray) = nothing
polarization!(ray::Ray, P) = nothing

# Placeholders for the future
intensity(ray::Ray) = nothing
intensity!(ray::Ray, I) = nothing


"""
    Ray(pos::Vector{P}, dir::Vector{D}, λ=1064) where {P, D}

Constructs an instance of `Ray` with common type `T`. Vector `dir` is normalized. Ray is initialized with `NoIntersection` and `λ` = 630 nm.\\
In general, `T<:AbstractFloat`! Inputs of type `Int` are promoted to `Float64`.
"""
function Ray(pos::Vector{P}, dir::Vector{D}, λ=1064) where {P<:Real, D<:Real}
    @assert norm3d(dir) != 0 "Illegal vector for direction"
    T = promote_type(P, D)
    # Catch integer type inputs
    if !(T<:AbstractFloat)
        T = Float64
    end
    return Ray{T}(T.(pos), normalize3d(T.(dir)), NoIntersection(T), Information(T.(λ)))
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

Temporary container struct to test ray tracing.
"""
mutable struct Beam{T}
    rays::Vector{Ray{T}}
end