"""
    Ray{T}

Mutable struct to store ray information. A `Ray` is described by ``\\vec{v}_{pos}+t\\cdot\\vec{v}_{dir}`` with ``t\\in[0,\\infty)``.

# Fields
- `id`: a UUID4 that uniquely identifies the `Ray`
- `pos`: a 3D-vector that describes the `Ray` origin
- `dir`: a normalized 3D-vector that describes the `Ray` direction
- `intersection`: refer to [`Intersection`](@ref)
- `parameters`: refer to [`Parameters`](@ref)
"""
mutable struct Ray{T} <: AbstractRay{T}
    id::UUID
    pos::Point3{T}
    dir::Point3{T}
    intersection::Nullable{Intersection{T}}
    parameters::Nullable{Parameters{T}}
end

function Base.length(ray::Ray{T}) where {T}
    if isnothing(intersection(ray))
        return Inf
    end
    return length(intersection(ray))
end

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

"""
    Ray(pos, dir, λ=1000e-9)

Constructs a `Ray` where:

- `pos`: is the `Ray` origin
- `dir`: is the `Ray` direction of propagation, normalized to unit length

Optionally, a wavelength `λ` can be specified.
"""
function Ray(pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ = 1000e-9) where {P <: Real, D <: Real}
    F = promote_type(P, D)
    return Ray{F}(uuid4(),
        Point3{F}(pos),
        normalize(Point3{F}(dir)),
        nothing,
        Parameters(λ))
end

"""
    line_point_distance3d(ray, point)

Returns value for the shortest distance between the `ray` (extended to ∞) and `point`.
"""
line_point_distance3d(ray::AbstractRay, point) = line_point_distance3d(position(ray),
    direction(ray),
    point)

"""
    angle3d(ray::AbstractRay, intersect::Intersection=intersection(ray))

Calculates the angle between a `ray` and its or some other `intersection`.
"""
function angle3d(ray::AbstractRay, intersect::Intersection = intersection(ray))
    return angle3d(direction(ray), normal3d(intersect))
end
