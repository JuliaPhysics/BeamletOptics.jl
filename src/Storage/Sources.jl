#=
Storage of sources: rays, beams, Gaussian beamlets and beam groups.

Only the untraced state is stored, i.e. the first ray of every beam. Traced rays, parents and
children are dropped, a loaded source starts without any intersection. Beam groups pack the
ray data column-wise into arrays, so that large groups become `.bin` assets.
=#

_check_source_eltype(x, ::Type{T}) where {T} =
    T === Float64 || throw(ArgumentError("only Float64 sources can be stored, got $(typeof(x))"))

#=
Rays
=#

function _ray_storage(r::AbstractRay{T}) where {T}
    _check_source_eltype(r, T)
    return Dict{String, Any}("pos" => encode_vec3(r.pos), "dir" => encode_vec3(r.dir),
        "λ" => Float64(r.λ), "n" => Float64(r.n))
end

to_storage(r::Ray, ctx) = _ray_storage(r)

function from_storage(::Type{Ray}, d, ctx)
    return Ray{Float64}(decode_vec3(d["pos"]), decode_vec3(d["dir"]), nothing,
        Float64(d["λ"]), Float64(d["n"]))
end

function to_storage(r::PolarizedRay, ctx)
    d = _ray_storage(r)
    d["E0"] = [encode_complex(e) for e in r.E0]
    return d
end

# The inner constructor normalizes `dir`, which may change its last bit: restore the stored value
function _polarized_ray(pos::Point3{Float64}, dir::Point3{Float64}, λ::Float64, n::Float64, E0::Point3{ComplexF64})
    ray = PolarizedRay{Float64}(pos, dir, nothing, λ, n, E0)
    ray.dir = dir
    return ray
end

function from_storage(::Type{PolarizedRay}, d, ctx)
    length(d["E0"]) == 3 || throw(ArgumentError("E0 of a PolarizedRay needs 3 entries"))
    E0 = Point3{ComplexF64}(decode_complex(d["E0"][1]), decode_complex(d["E0"][2]), decode_complex(d["E0"][3]))
    return _polarized_ray(decode_vec3(d["pos"]), decode_vec3(d["dir"]),
        Float64(d["λ"]), Float64(d["n"]), E0)
end

register_storage_type!(Ray, "Ray")
register_storage_type!(PolarizedRay, "PolarizedRay")

#=
Packed rays of groups: one column per ray
=#

const _PACKED_RAY_TYPES = Dict{String, Type}("Ray" => Ray{Float64}, "PolarizedRay" => PolarizedRay{Float64})

_ray_eltype(::Type{<:AbstractRay{T}}) where {T} = T

function _pack_rays(::Type{R}, rs::AbstractVector, ctx; stem::AbstractString) where {R <: AbstractRay}
    _check_source_eltype(rs, _ray_eltype(R))
    R <: Union{Ray, PolarizedRay} || throw(ArgumentError("rays of type $R cannot be stored"))
    polarized = R <: PolarizedRay
    N = length(rs)
    pos = Matrix{Float64}(undef, 3, N)
    dir = Matrix{Float64}(undef, 3, N)
    λ = Vector{Float64}(undef, N)
    n = Vector{Float64}(undef, N)
    E0 = Matrix{ComplexF64}(undef, 3, polarized ? N : 0)
    for (i, r) in enumerate(rs)
        pos[:, i] = r.pos
        dir[:, i] = r.dir
        λ[i] = r.λ
        n[i] = r.n
        polarized && (E0[:, i] = r.E0)
    end
    d = Dict{String, Any}("ray_type" => polarized ? "PolarizedRay" : "Ray", "count" => N,
        "pos" => encode_array(pos, ctx; stem = "$stem-pos"),
        "dir" => encode_array(dir, ctx; stem = "$stem-dir"),
        "λ" => encode_array(λ, ctx; stem = "$stem-lambda"),
        "n" => encode_array(n, ctx; stem = "$stem-n"))
    polarized && (d["E0"] = encode_array(E0, ctx; stem = "$stem-E0"))
    return d
end

function _unpack_rays(d, ctx)
    R = get(_PACKED_RAY_TYPES, d["ray_type"]) do
        throw(ArgumentError("unknown ray type \"$(d["ray_type"])\""))
    end
    N = Int(d["count"])
    pos = decode_array(d["pos"], ctx)
    dir = decode_array(d["dir"], ctx)
    λ = decode_array(d["λ"], ctx)
    n = decode_array(d["n"], ctx)
    mismatch() = throw(ArgumentError("packed ray data does not match count $N"))
    (size(pos) == size(dir) == (3, N) && length(λ) == length(n) == N) || mismatch()
    col(A, i) = Point3{eltype(A)}(A[1, i], A[2, i], A[3, i])
    if R <: PolarizedRay
        E0 = decode_array(d["E0"], ctx)
        size(E0) == (3, N) || mismatch()
        return PolarizedRay{Float64}[_polarized_ray(col(pos, i), col(dir, i), λ[i], n[i], col(E0, i)) for i in 1:N]
    end
    return Ray{Float64}[Ray{Float64}(col(pos, i), col(dir, i), nothing, λ[i], n[i]) for i in 1:N]
end

#=
Beams
=#

to_storage(b::Beam, ctx) = Dict{String, Any}("ray" => encode(first(rays(b)), ctx))

from_storage(::Type{Beam}, d, ctx) = Beam(decode(d["ray"], ctx))

register_storage_type!(Beam, "Beam")

function to_storage(g::GaussianBeamlet{T}, ctx) where {T}
    _check_source_eltype(g, T)
    return Dict{String, Any}(
        "chief" => encode(first(rays(g.chief)), ctx),
        "waist" => encode(first(rays(g.waist)), ctx),
        "divergence" => encode(first(rays(g.divergence)), ctx),
        "λ" => Float64(g.λ), "w0" => Float64(g.w0), "E0" => encode_complex(g.E0))
end

function from_storage(::Type{GaussianBeamlet}, d, ctx)
    beam(key) = Beam(decode(d[key], ctx)::Ray{Float64})
    return GaussianBeamlet(beam("chief"), beam("waist"), beam("divergence"),
        Float64(d["λ"]), Float64(d["w0"]), decode_complex(d["E0"]))
end

register_storage_type!(GaussianBeamlet, "GaussianBeamlet")

const _ASTIGMATIC_AUX = (:wxp, :wxm, :wyp, :wym, :dxp, :dxm, :dyp, :dym)

function to_storage(agb::AstigmaticGaussianBeamlet{T}, ctx) where {T}
    _check_source_eltype(agb, T)
    d = Dict{String, Any}("c" => encode(first(rays(agb.c)), ctx))
    for key in _ASTIGMATIC_AUX
        d[string(key)] = encode(first(rays(getfield(agb, key))), ctx)
    end
    return d
end

function from_storage(::Type{AstigmaticGaussianBeamlet}, d, ctx)
    c = Beam(decode(d["c"], ctx)::PolarizedRay{Float64})
    aux = (Beam(decode(d[string(key)], ctx)::Ray{Float64}) for key in _ASTIGMATIC_AUX)
    return AstigmaticGaussianBeamlet(c, aux...)
end

register_storage_type!(AstigmaticGaussianBeamlet, "AstigmaticGaussianBeamlet")

#=
Beam groups. `UniformDiscSource` returns a `CollimatedSource` and needs no tag of its own.
=#

_pack_beams(bs::Vector{Beam{T, R}}, ctx) where {T, R} =
    _pack_rays(R, [first(rays(b)) for b in bs], ctx; stem = "beams")

_unpack_beams(d, ctx) = [Beam(r) for r in _unpack_rays(d, ctx)]

function to_storage(ps::PointSource{T}, ctx) where {T}
    _check_source_eltype(ps, T)
    return Dict{String, Any}("beams" => _pack_beams(ps.beams, ctx), "NA" => Float64(ps.NA))
end

from_storage(::Type{PointSource}, d, ctx) = PointSource(_unpack_beams(d["beams"], ctx), Float64(d["NA"]))

register_storage_type!(PointSource, "PointSource")

function to_storage(cs::CollimatedSource{T}, ctx) where {T}
    _check_source_eltype(cs, T)
    return Dict{String, Any}("beams" => _pack_beams(cs.beams, ctx), "diameter" => Float64(cs.diameter))
end

from_storage(::Type{CollimatedSource}, d, ctx) =
    CollimatedSource(_unpack_beams(d["beams"], ctx), Float64(d["diameter"]))

register_storage_type!(CollimatedSource, "CollimatedSource")

# Chief rays in one packed table, the 8 auxiliary rays of each beamlet in a second one (beamlet-major)
function to_storage(bg::AstigmaticBeamGroup{T}, ctx) where {T}
    _check_source_eltype(bg, T)
    chief = [first(rays(b.c)) for b in bg.beams]
    aux = [first(rays(getfield(b, key))) for b in bg.beams for key in _ASTIGMATIC_AUX]
    return Dict{String, Any}(
        "chief" => _pack_rays(PolarizedRay{T}, chief, ctx; stem = "chief"),
        "aux" => _pack_rays(Ray{T}, aux, ctx; stem = "aux"))
end

function from_storage(::Type{AstigmaticBeamGroup}, d, ctx)
    chief = _unpack_rays(d["chief"], ctx)
    aux = _unpack_rays(d["aux"], ctx)
    (eltype(chief) === PolarizedRay{Float64} && eltype(aux) === Ray{Float64} &&
     length(aux) == 8 * length(chief)) ||
        throw(ArgumentError("AstigmaticBeamGroup needs polarized chief rays and 8 auxiliary rays per beamlet"))
    beamlets = AstigmaticGaussianBeamlet{Float64}[
        AstigmaticGaussianBeamlet(Beam(c), (Beam(aux[8(i - 1) + k]) for k in 1:8)...)
        for (i, c) in enumerate(chief)]
    return AstigmaticBeamGroup(beamlets)
end

register_storage_type!(AstigmaticBeamGroup, "AstigmaticBeamGroup")
