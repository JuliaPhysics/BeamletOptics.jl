"""
    Interaction{T}

Describes how a `ray` and an `object` interact with each other.
Returns the new parameters for the construction of the next ray.
"""
struct Interaction{T}
    pos::Vector{T}
    dir::Vector{T}
    parameters::Parameters{T}
    hint::Nullable{UUID}
end

position(interaction::Interaction) = interaction.pos

direction(interaction::Interaction) = interaction.dir

parameters(interaction::Interaction) = interaction.parameters

hint(interaction::Interaction) = interaction.hint

"""
    AbstractReflectiveOptic <: AbstractObject

A generic type to represent `AbstractObject`s which reflect incoming rays. The main function of `interact3d` should be akin to [`reflection3d`](@ref).
"""
abstract type AbstractReflectiveOptic <: AbstractObject end

function interact3d(::AbstractReflectiveOptic, ray::AbstractRay{R}) where {R<:Real}
    normal = intersection(ray).n
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    return Interaction{R}(npos, ndir, ray.parameters, nothing)
end

"""
    Mirror{S <: AbstractShape} <: AbstractReflectiveOptic

Concrete implementation of a perfect mirror with arbitrary shape.
"""
struct Mirror{S<:AbstractShape} <: AbstractReflectiveOptic
    id::UUID
    shape::S
end

"""
    AbstractRefractiveOptic <: AbstractObject

A generic type to represent `AbstractObject`s which refract incoming rays. The main function of `interact3d` should be akin to [`refraction3d`](@ref).

# Implementation reqs.
Subtypes of `AbstractRefractiveOptic` should implement all supertype reqs as well as:

# Fields
- `n::Function`: a function which returns the refractive index for a wavelength λ
"""
abstract type AbstractRefractiveOptic <: AbstractObject end

refractive_index(object::AbstractRefractiveOptic) = object.n

function interact3d(object::AbstractRefractiveOptic, ray::Ray{R}) where {R}
    # Check dir. of ray and surface normal
    normal = intersection(ray).n
    λ = wavelength(ray)
    if dot(direction(ray), normal) < 0
        # "Outside prism"
        n1 = 1.0
        n2 = refractive_index(object)(λ)
    else
        # "Inside prism"
        n1 = refractive_index(object)(λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    # Hint is the current object ID
    return Interaction{R}(npos, ndir, Parameters(λ, n2), id(object))
end

struct Lens{S<:AbstractShape,T<:Function} <: AbstractRefractiveOptic
    id::UUID
    shape::S
    n::T
end

SphericalLens(r1, r2, l, d, n) = SphericalLens(r1, r2, l, d, λ->n)

"""
    SphericalLens(r1, r2, l=0, d=1SCDI.inch, n=λ->1.5)

Creates a spherical lens based on:

- `r1`: front radius
- `r2`: back radius
- `l`: length of cylindrical section, default is 0
- `d`: lens diameter, default is one inch
- `n`: refractive index as a function of λ

To differentiate surface types:

- `r > 0`: convex
- `r < 0`: concave
- `r = Inf`: planar
"""
function SphericalLens(r1::Real, r2::Real, l::Real=0, d::Real=1SCDI.inch, n::Function=λ->1.5)
    if !isfinite(r1) || isnan(r1)
        error("r1 can not be Inf")
    end
    # Determine lens SDF
    shape::SCDI.Nullable{SCDI.AbstractSphericalLensSDF} = nothing
    # Test for thin lens
    if iszero(l)
        shape = SCDI.ThinLensSDF(r1, r2, d)
        # goto to avoid overwrite
        @goto lens_creator
    end
    # Test for plano lens
    if isinf(r2)
        if r1 > 0
            shape = SCDI.PlanoConvexLensSDF(r1, l, d)
        elseif r1 < 0
            shape = SCDI.PlanoConcaveLensSDF(abs(r1), l, d)
        end
        @goto lens_creator
    end
    # Test for bi-convex/concave or meniscus
    if r1 > 0 && r2 > 0
        shape = SCDI.BiConvexLensSDF(r1, r2, l, d)
    end
    if r1 < 0 && r2 < 0
        shape = SCDI.BiConcaveLensSDF(abs(r1), abs(r2), l, d)
    end
    if r1 > 0 && r2 < 0
        shape = SCDI.ConvexConcaveLensSDF(r1, abs(r2), l, d)
    end
    # Create lens
    @label lens_creator
    if isnothing(shape)
        error("Could not find suitable lens SDF")
    end
    return SCDI.Lens(uuid4(), shape, n)
end

ThinLens(R1::Real, R2::Real, d::Real, n::Real) = ThinLens(R1, R2, d, x -> n)
function ThinLens(R1::Real, R2::Real, d::Real, n::Function)
    shape = ThinLensSDF(R1, R2, d)
    return Lens(uuid4(), shape, n)
end

struct Prism{S<:AbstractShape,T<:Function} <: AbstractRefractiveOptic
    id::UUID
    shape::S
    n::T
end

#=
Implements photodetector, efield calculation after solve_system! with photodetection!
=#
abstract type AbstractDetector <: AbstractObject end

mutable struct Photodetector{S<:AbstractShape,T} <: AbstractDetector
    const id::UUID
    const shape::S
    rays::Nullable{Vector{UUID}}
    x::LinRange{T,Int64}
    y::LinRange{T,Int64}
    field::Matrix{Complex{T}}
end

function Photodetector(scale::T, n::Int) where {T}
    sz = 0.5
    vertices = [
        sz 0 sz
        sz 0 -sz
        -sz 0 -sz
        -sz 0 sz
    ]
    faces = [
        1 2 4
        2 3 4
    ]
    x = y = LinRange(-sz, sz, n) * scale
    field = zeros(Complex{T}, n, n)
    shape = Mesh{T}(uuid4(), vertices .* scale, faces, Matrix{T}(I, 3, 3), T.([0, 0, 0]), scale)
    return Photodetector{typeof(shape),T}(uuid4(), shape, nothing, x, y, field)
end

function interact3d(pd::Photodetector, ray::Ray)
    if isnothing(pd.rays)
        pd.rays = [id(ray)]
    else
        push!(pd.rays, id(ray))
    end
    return nothing
end

function photodetection!(pd::Photodetector, gauss::GaussianBeamlet)
    # Test if beamlet hits pd
    is_hit_by = false
    for ray_id in pd.rays
        if !isparentbeam(gauss, ray_id)
            continue
        end
        if ray_id == last(gauss.chief.rays).id
            is_hit_by = true
            break
        end
    end
    # Exit if not hit by gauss
    if !is_hit_by
        return nothing
    end
    # Add efield contribution to pd.field
    ray = gauss.chief.rays[end]
    l = length(gauss.chief)
    T = transpose(orientation(shape(pd)))
    P = ray.pos + ray.dir * ray.intersection.t
    shape_pos = position(shape(pd))
    Threads.@threads for (j, y) in collect(enumerate(pd.y))     # row column major order?
        # thread-local storage
        p = Vector{Float64}(undef, 3)
        # the second entry does not depend on x
        p[2] = shape_pos[2]
        # reusable line point buffers
        lp_a = Vector{Float64}(undef, 3)
        lp_b = Vector{Float64}(undef, 3)
        # reusable isinfrontof buffer
        los = Vector{Float64}(undef, 3)
        # reusable gauss/point buffer
        g_b = Vector{Float64}(undef, 3)
        p_b = Vector{Float64}(undef, 3)
        @inbounds for (i, x) in enumerate(pd.x)
            # Transform point p on PD into world coords
            p[1] = T[1, 1] * x + T[1, 3] * y + shape_pos[1]
            p[3] = T[3, 1] * x + T[3, 3] * y + shape_pos[3]
            r = line_point_distance3d(ray, p, lp_a, lp_b)
            c = sqrt(sum(x -> (x[1] - x[2])^2, zip(P, p)))
            z = sqrt(abs(c^2 - r^2)) # abs to protect against small neg. values
            # Correct sign of z
            if isinfrontof(p, P, ray.dir, los)
                z = -z
            end
            pd.field[i, j] += electric_field(gauss, r, l + z, g_b, p_b)
        end
    end
    return nothing
end

intensity(pd::Photodetector) = intensity.(pd.field)

"""
    optical_power(pd::Photodetector)

Calculates the total optical power on `pd` in [W] by integration over the local intensity.
"""
optical_power(pd::Photodetector) = trapz((pd.x, pd.y), intensity(pd))

reset_photodetector!(pd::Photodetector{S,T}) where {S,T} = (pd.field .= zero(Complex{T}))

"""
    photodetector_resolution!(pd::Photodetector, n::Int)

Sets the resolution of `pd` to `n` × `n`. Note that this resets the current `pd.field`.
"""
function photodetector_resolution!(pd::Photodetector{S,T}, n::Int) where {S,T}
    pd.x = LinRange(pd.x.start, pd.x.stop, n)
    pd.y = LinRange(pd.y.start, pd.y.stop, n)
    pd.field = zeros(Complex{T}, n, n)
    return nothing
end