"""
    AbstractBeamGroup

FIXME
"""
abstract type AbstractBeamGroup{T <: Real, R <: AbstractRay{T}} end

beams(bg::AbstractBeamGroup) = bg.beams

position(bg::AbstractBeamGroup) = position(first(rays(first(beams(bg)))))
direction(bg::AbstractBeamGroup) = direction(first(rays(first(beams(bg)))))

function solve_system!(system::AbstractSystem, bg::AbstractBeamGroup; kwargs...)
    for beam in beams(bg)
        solve_system!(system, beam; kwargs...)
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", bg::AbstractBeamGroup)
    println(io, "Subtype of AbstractBeamGroup")
    println(io, "   # of beams: $(length(beams(bg)))")
    return nothing
end

struct PointSource{T, R<:AbstractRay{T}} <: AbstractBeamGroup{T,R}
    NA::T
    beams::Vector{Beam{T, R}}
end

"""
    PointSource

FIXME
"""
function PointSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ::H,
        λ::L;
        num_rays::Int=1000,
        num_rings::Int=10
    ) where {P, D, H, L}
    T = promote_type(P, D, H, L)
    if num_rays < num_rings
        @warn "# of rays should be > than # of rings"
    end
    if θ ≥ pi
        throw(ErrorException("Point source openign half-angle θ must be ≤ π"))
    end
    # define basis vectors
    dir = normalize(dir)
    b1 = normal3d(dir)
    b2 = normal3d(dir, b1)
    θ_NA = LinRange(0, θ, num_rings)
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # calculate total accumulated circumference of all rings
    ndirs = [rotate3d(b2, step(θ_NA)*i) * dir for i in eachindex(θ_NA[2:end])]
    circm = norm.(ndirs .- dot.(ndirs, Ref(dir)) .* Ref(dir)) .* Ref(2π)
    total = sum(circm)
    ds = total / num_rays
    # calculate number of rays per ring
    n_rays = round.(Int, circm / ds)
    # correct n_rays to match num_rays
    n_rays[end] += (num_rays - sum(n_rays))
    for (i, ndir) in enumerate(ndirs)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir, dphi)
        cdir = ndir
        for _ = 1:numEl
            push!(beams, Beam(pos, cdir, λ))
            # rotate vector (not-thread safe!)
            cdir = RotMat * cdir
        end
    end
    NA = numerical_aperture(θ)
    return PointSource(NA, beams)
end

struct CollimatedSource{T, R<:AbstractRay{T}} <: AbstractBeamGroup{T,R}
    beams::Vector{Beam{T, R}}
end

"""
    CollimatedSource

FIXME
"""
function CollimatedSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        diameter::D2,
        λ::L;
        num_rings::Int=20,
        num_rays::Int=1000
    ) where {P, D1, D2, L}
    T = promote_type(P, D1, D2, L)
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # setup concentric beam ring radii
    b1 = normal3d(dir)
    r_max = diameter/2
    radii = LinRange(0, r_max, num_rings)[2:end]
    # calculate total accumulated circumference of all rings
    circm = radii*2π
    total = sum(circm)
    ds = total / num_rays
    # calculate number of rays per ring
    n_rays = round.(Int, circm / ds)
    # correct n_rays to match num_rays
    n_rays[end] += (num_rays - sum(n_rays))
    # Generate beam rings
    for (i, r) in enumerate(radii)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir, dphi)
        helper = b1 * r
        for _ = 1:numEl
            push!(beams, Beam(pos + helper, dir, λ))
            helper = RotMat * helper
        end
    end
    return CollimatedSource(beams)
end

export PointSource, CollimatedSource