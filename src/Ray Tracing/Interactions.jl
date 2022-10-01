struct Interaction{T}
    pos::Union{Vector{T}, Missing}
    dir::Union{Vector{T}, Missing}
    information::Information{T}
end

# Optimized for pointer look-up
const _NoInteractionF64 = Interaction{Float64}(missing, missing, NoInformation(Float64))
const _NoInteractionF32 = Interaction{Float32}(missing, missing, NoInformation(Float32))

NoInteraction(::Type{Float64}) = _NoInteractionF64
NoInteraction(::Type{Float32}) = _NoInteractionF32

struct Mirror{T} <: AbstractMesh{T}
    mesh::Mesh{T}
end

function interact3d(mirror::Mirror{M}, ray::AbstractRay{R}) where {M, R}
    T = promote_type(M, R)
    normal = ray.intersection.n
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    return Interaction{T}(npos, ndir, ray.information)
end

struct Prism{T} <: AbstractMesh{T}
    mesh::Mesh{T}
    ref_index::Function
end

refractive_index(prism::Prism) = prism.ref_index
refractive_index!(prism, ref_index) = nothing

function interact3d(prism::Prism{P}, ray::Ray{R}) where {P, R}
    T = promote_type(P, R)
    # Check dir. of ray and surface normal
    normal = ray.intersection.n
    λ = wavelength(ray)
    if dot(direction(ray), normal) < 0
        # "Outside prism"
        n1 = 1.0
        n2 = refractive_index(prism)(λ)
    else
        # "Inside prism"
        n1 = refractive_index(prism)(λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    return Interaction{T}(npos, ndir, Information(λ, n2))
end

function interact3d(lens::SingletLens{G,V,M,W,F}, ray::Ray{R}) where {G,V,M,W,F,R}
    T = promote_type(W, R)
    # Check dir. of ray and surface normal
    normal = ray.intersection.n
    λ = wavelength(ray)
    if dot(direction(ray), normal) < 0
        # "Outside lens"
        n1 = refractive_index(ray)
        n2 = refractive_index(lens)(λ)
    else
        # "Inside lens"
        n1 = refractive_index(lens)(λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    return Interaction{T}(npos, ndir, Information(λ, n2))
end

refractive_index(lens::SingletLens) = lens.ref_index
