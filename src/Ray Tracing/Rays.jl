"""
    Intersection{T}

Stores some data resulting from ray tracing a `System`. This information can be used, i.e. for retracing.

# Fields:
- `t`: length of the ray parametrization in [m]
- `n`: normal vector at the point of intersection
- `id`: index of the intersected object in the `System` object-vector
"""
mutable struct Intersection{T} <: AbstractEntity
    t::T
    n::Vector{T}
    id::Nullable{UUID}
end


Base.length(intersection::Intersection) = intersection.t
length!(::Intersection, t) = nothing

normal3d(in::Intersection) = in.n

"""
    Parameters{T}

Stores the optical parameters that are relevant for the propagation of a ray through an optical system. See also `Ray{T}`.

# Fields:
- `λ`: wavelength in [nm]
- `n`: refractive index along the beam path
- `I`: intensity value [**currently a placeholder**]
- `P`: polarization vector (either Stokes or Jones formalism, tbd.) [**currently a placeholder**]
"""
mutable struct Parameters{T}
    λ::T
    n::T
    I::T
    P::Vector{Complex{T}}
end

"Prototype constructor for `Parameters{T}`"
function Parameters(λ::T, n=1, I=1, P=[0, 0]) where {T}
    return Parameters{T}(λ, n, I, P)
end

"""
    Ray{T}

Mutable struct to store ray information. A `Ray` is described by ``\\vec{v}_{pos}+t\\cdot\\vec{v}_{dir}`` with ``t\\in[0,\\infty)``.

# Fields
- `pos`: a 3D-vector that describes the `Ray` origin
- `dir`: a normalized 3D-vector that describes the `Ray` direction
- `intersection`: refer to `Intersection{T}`.
- `parameters`: refer to `Parameters{T}`
"""
mutable struct Ray{T} <: AbstractRay{T}
    id::UUID
    pos::Vector{T}
    dir::Vector{T}
    intersection::Nullable{Intersection{T}}
    parameters::Nullable{Parameters{T}}
end

Base.length(ray::Ray{T}) where T = ray.intersection.t
length!(::Ray, len) = nothing

intersection(ray::Ray) = ray.intersection
intersection!(ray::Ray, intersection) = (ray.intersection = intersection)

parameters(ray::Ray) = ray.parameters
parameters!(ray::Ray, parameters) = (ray.parameters = parameters)

wavelength(ray::Ray) = ray.parameters.λ
wavelength!(ray::Ray, λ) = (ray.parameters.λ = λ)

refractive_index(ray::Ray) = ray.parameters.n
refractive_index!(ray::Ray, n) = (ray.parameters.n = n)

# Placeholders for the future
polarization(::Ray) = nothing
polarization!(::Ray, P) = nothing

# Placeholders for the future
intensity(ray::Ray) = ray.parameters.I
intensity!(ray::Ray, I) = (ray.parameters.I = I)

function Ray(pos::Vector{T}, dir::Vector{T}, λ=1064) where {T<:AbstractFloat}
    return Ray{T}(uuid4(), pos, normalize3d(dir), nothing, Parameters(T(λ)))
end

"""
    Ray(pos::Vector{P}, dir::Vector{D}, λ=1064) where {P, D}

Constructs an instance of `Ray` with common type `T`. Vector `dir` is normalized. Ray is initialized with no `intersection` and `λ` = 1064 nm.\\
In general, `T<:AbstractFloat`! Inputs of type `Int` are promoted to `Float64`.
"""
function Ray(pos::Vector{P}, dir::Vector{D}, λ=1064) where {P<:Real,D<:Real}
    F = float(promote_type(P, D))

    return Ray(convert(Vector{F}, pos), convert(Vector{F}, dir), λ)
end

"""
    line_point_distance3d(ray, point)

Returns value for the shortest distance between the `ray` (extended to ∞) and `point`.
"""
line_point_distance3d(ray::AbstractRay, point, buf_a=zeros(length(point)), buf_b=zeros(length(point))) = line_point_distance3d(position(ray), direction(ray), point, buf_a, buf_b)

"""
    angle3d(ray::AbstractRay, intersect::Intersection=intersection(ray))

Calculates the angle between a `ray` and its or some other `intersection`.
"""
function angle3d(ray::AbstractRay, intersect::Intersection=intersection(ray))
    return angle3d(direction(ray), normal3d(intersect))
end

"""
    Beam

Temporary container struct to test ray tracing.
"""
mutable struct Beam{T}
    rays::Vector{Ray{T}}
end

Beam(ray::Ray{T}) where T = Beam{T}([ray])

"""
    length(beam::Beam)

Calculate the length of a beam up to the point of the last intersection.
"""
function Base.length(beam::Beam{T}) where T
    l = zero(T)
    for ray in beam.rays
        if isnothing(intersection(ray))
            break
        end
        l += length(ray)
    end
    return l
end

"""
    point_on_beam(beam::Beam, t::Real)

Function to find a point given a specific distance `t` along the beam. Return the ray `index` aswell.
For negative distances, assume first ray backwards.
"""
point_on_beam(beam::Beam{T}, t::Real) where T = point_on_beam!(Vector{T}(undef, 3), beam, t)

function point_on_beam!(point::AbstractVector, beam::Beam, t::Real)
    # Initialize counter to track cumulative length
    temp = zero(t)
    numEl = length(beam.rays)
    for (index, ray) in enumerate(beam.rays)
        # Catch final ray
        if index == numEl
            break
        end
        temp += length(ray)
        # If the specified distance `t` is less than the cumulative length,
        # calculate the local ray length `b` and find the point along the ray
        if t < temp
            b = temp - t
            @. point = $position(ray) + ($length(ray) - b) * $direction(ray)
            return point, index
        end
    end
    # If no solution at this point assume final ray with infinite length
    ray = beam.rays[end]
    b = t - temp
    @. point .= $position(ray) + b * $direction(ray)
    return point, length(beam.rays)
end

"""
    isparaxial(system, beam, threshold=π/4)

Tests the angle between the `beam` direction and surface normal at each intersection.
Mainly intended as a check for [`GaussianBeamlet`](@ref).
"""
function isparaxial(system::AbstractSystem, beam::Beam, threshold::Real=π/4)
    # Test if refractive elements are hit with angle larger than threshold
    for ray in beam.rays
        if isnothing(intersection(ray))
            break
        end
        target = object(system, id(intersection(ray)))
        # Test if refractive element
        if !isa(target, AbstractRefractiveOptic)
            continue
        end
        # Test angle between ray and its intersection
        angle = angle3d(ray)
        if angle > π/2 # flip sector
            angle = π - angle
        end
        if angle > threshold # rad
            return false
        end
    end
    return true
end

function isparentbeam(beam::Beam, ray_id::UUID)
    for ray in beam.rays
        if id(ray) == ray_id
            return true
        end
    end
    return false
end

"""
    isparentbeam(beam, ray)

Tests if the given `beam` contains the `ray` as a part of its solution.
"""
isparentbeam(beam::Beam, ray::AbstractRay) = isparentbeam(beam, ray.id)

function Base.show(io::IO, ::MIME"text/plain", beam::Beam)
    for (i, ray) in enumerate(beam.rays)
        println(io, "Ray $i:")
        if isnothing(intersection(ray))
            println(io, "    No intersection")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
        else
            println(io, "    Intersects with object #$(intersection(ray).id)")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
            println(io, "    End.: $(position(ray) .+ length(ray) .* direction(ray))")
        end
    end
    return nothing
end