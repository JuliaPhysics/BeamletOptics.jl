# ==============================================================================
# FourierBeamletPropagator.jl
#
# Part of the Fourier-optics detector architecture in BeamletOptics.jl.
# Provides a lightweight, high-performance Fourier beamlet propagator that
# transforms incident rays and Gaussian beamlets into discrete Fourier k-points
# (FourierKPoint) for frequency-space propagation and evaluation.
# ==============================================================================

"""
    FourierKPoint{T <: Real} <: AbstractDetectorHit

A discrete frequency-space point representing an electromagnetic plane wave component.

Parameters:
- `k4`: Four-wavevector `[kx, ky, kz, omega]` of type `SVector{4, T}` with spatial wavevector components in rad/m and angular frequency in rad/s.
- `E`: Complex 3D electric field vector `[Ex, Ey, Ez]` of type `SVector{3, Complex{T}}`.
- `weight`: Quadrature / spectral weight of type `T`.

Contracts & Invariants:
- Angular optical frequency `k4[4]` (or `omega`) is non-negative for propagating plane waves.
- Backward compatibility is fully preserved via property accessors `p.k` and `p.omega`.
"""
struct FourierKPoint{T <: Real} <: AbstractDetectorHit
    k4::SVector{4, T} # [kx, ky, kz, omega]
    E::SVector{3, Complex{T}}
    weight::T

    @inline function FourierKPoint{T}(k4::SVector{4, T}, E::SVector{3, Complex{T}}, weight::T = one(T)) where {T <: Real}
        new{T}(k4, E, weight)
    end
end

@inline FourierKPoint(k4::SVector{4, T}, E::SVector{3, Complex{T}}, weight::T = one(T)) where {T <: Real} =
    FourierKPoint{T}(k4, E, weight)

"""
    Base.getproperty(p::FourierKPoint{T}, s::Symbol) where {T}

Backward-compatible property accessor providing access to `:k` and `:omega`.

Parameters:
- `p`: Fourier k-point instance of type `FourierKPoint{T}`.
- `s`: Property name symbol.

Returns:
- 3-vector wavevector if `s === :k`, angular frequency scalar if `s === :omega`, or direct field value.

Contracts & Invariants:
- Inlined and allocation-free (@allocated == 0).
"""
@inline Base.getproperty(p::FourierKPoint{T}, s::Symbol) where {T} =
    s === :k ? SVector{3, T}(getfield(p, :k4)[1], getfield(p, :k4)[2], getfield(p, :k4)[3]) :
    s === :omega ? getfield(p, :k4)[4] : getfield(p, s)

@inline Base.propertynames(::FourierKPoint, private::Bool = false) = (:k4, :E, :weight, :k, :omega)

"""
    wavenumber(p::FourierKPoint{T}) where {T}

Returns the magnitude of the spatial wavevector in rad/m.

Parameters:
- `p`: Fourier k-point instance of type `FourierKPoint{T}`.

Returns:
- Scalar wavevector magnitude of type `T` in rad/m.
"""
@inline wavenumber(p::FourierKPoint{T}) where {T} = norm(p.k)

"""
    wavelength(p::FourierKPoint{T}) where {T}

Returns the optical wavelength in meters.

Parameters:
- `p`: Fourier k-point instance of type `FourierKPoint{T}`.

Returns:
- Scalar wavelength of type `T` in meters (defaults to `one(T)` if wavenumber is zero).
"""
@inline wavelength(p::FourierKPoint{T}) where {T} = (k_mag = norm(p.k); k_mag > 0 ? (T(2π) / k_mag) : one(T))

"""
    direction(p::FourierKPoint{T}) where {T}

Returns the normalized wave propagation direction unit 3-vector.

Parameters:
- `p`: Fourier k-point instance of type `FourierKPoint{T}`.

Returns:
- Unit direction vector of type `SVector{3, T}`.
"""
@inline direction(p::FourierKPoint{T}) where {T} = (k_mag = norm(p.k); k_mag > 0 ? (p.k / k_mag) : SVector{3, T}(zero(T), one(T), zero(T)))

"""
    position(p::FourierKPoint{T}) where {T}

Returns the spatial reference position 3-vector in meters.

Parameters:
- `p`: Fourier k-point instance of type `FourierKPoint{T}`.

Returns:
- Spatial reference position of type `SVector{3, T}` in meters.
"""
@inline position(p::FourierKPoint{T}) where {T} = SVector{3, T}(zero(T), zero(T), zero(T))

"""
    polarization(p::FourierKPoint)

Returns the complex 3D electric field polarization vector.

Parameters:
- `p`: Fourier k-point instance of type `FourierKPoint`.

Returns:
- Complex electric field 3-vector in V/m.
"""
@inline polarization(p::FourierKPoint) = p.E

"""
    FourierKPoint{T}(k4::AbstractVector{<:Real}, E::AbstractVector{<:Number}, weight::Real = one(T)) where {T <: Real}

Constructs a `FourierKPoint{T}` from 4-vector `k4`, complex electric field `E`, and optional weight.
"""
@inline function FourierKPoint{T}(k4::AbstractVector{<:Real}, E::AbstractVector{<:Number}, weight::Real = one(T)) where {T <: Real}
    return FourierKPoint{T}(
        SVector{4, T}(T(k4[1]), T(k4[2]), T(k4[3]), T(k4[4])),
        SVector{3, Complex{T}}(Complex{T}(E[1]), Complex{T}(E[2]), Complex{T}(E[3])),
        T(weight)
    )
end

"""
    FourierKPoint(k4::AbstractVector{<:Real}, E::AbstractVector{<:Number}, weight::Real = 1.0)

Constructs a `FourierKPoint` from 4-vector `k4`, complex electric field `E`, and optional weight with type promotion.
"""
@inline function FourierKPoint(k4::AbstractVector{<:Real}, E::AbstractVector{<:Number}, weight::Real = 1.0)
    T = promote_type(eltype(k4), real(eltype(E)), typeof(float(weight)))
    return FourierKPoint{T}(k4, E, weight)
end

"""
    FourierKPoint{T}(k::AbstractVector{<:Real}, omega::Real, E::AbstractVector{<:Number}, weight::Real = one(T)) where {T <: Real}

Constructs a `FourierKPoint{T}` from 3-vector wavevector `k`, angular frequency `omega`, complex electric field `E`, and optional weight.

Parameters:
- `k`: Spatial wavevector of type `AbstractVector{<:Real}` in rad/m.
- `omega`: Angular optical frequency of type `Real` in rad/s.
- `E`: Complex 3D electric field vector of type `AbstractVector{<:Number}`.
- `weight`: Spectral weight of type `Real`, defaults to `one(T)`.

Returns:
- Newly constructed `FourierKPoint{T}` with 4-vector `k4`.

Contracts & Invariants:
- Preserves type stability and target precision `T`.
"""
@inline function FourierKPoint{T}(k::AbstractVector{<:Real}, omega::Real, E::AbstractVector{<:Number}, weight::Real = one(T)) where {T <: Real}
    return FourierKPoint{T}(
        SVector{4, T}(T(k[1]), T(k[2]), T(k[3]), T(omega)),
        SVector{3, Complex{T}}(Complex{T}(E[1]), Complex{T}(E[2]), Complex{T}(E[3])),
        T(weight)
    )
end

"""
    FourierKPoint(k::AbstractVector{<:Real}, omega::Real, E::AbstractVector{<:Number}, weight::Real = 1.0)

Constructs a `FourierKPoint` by promoting numerical argument types to a common real precision type `T`.

Parameters:
- `k`: Spatial wavevector of type `AbstractVector{<:Real}` in rad/m.
- `omega`: Angular optical frequency of type `Real` in rad/s.
- `E`: Complex 3D electric field vector of type `AbstractVector{<:Number}`.
- `weight`: Spectral weight of type `Real`, defaults to 1.0.

Returns:
- `FourierKPoint{T}` with promoted real type `T`.

Contracts & Invariants:
- Promotes all types consistently.
"""
@inline function FourierKPoint(k::AbstractVector{<:Real}, omega::Real, E::AbstractVector{<:Number}, weight::Real = 1.0)
    T = promote_type(eltype(k), typeof(float(omega)), real(eltype(E)), typeof(float(weight)))
    return FourierKPoint{T}(k, omega, E, weight)
end

"""
    AbstractFourierDetector{T} <: AbstractDetector{T}

Abstract supertype for all Fourier-optics detectors based on frequency-space interpolation
and multidimensional fast Fourier transforms.
"""
abstract type AbstractFourierDetector{T} <: AbstractDetector{T} end


"""
    FourierBeamletPropagator{T <: Real, S <: AbstractShape{T}} <: AbstractFourierDetector{T}

A lean Fourier detector and propagator that captures ray and beamlet hits and
evaluates their discrete angular spectrum `Vector{FourierKPoint{T}}`.

# Fields
- `shape`: Geometric boundary shape.
- `hits`: Captured detector hit records.
- `stop`: If `true`, absorbs incoming beams; if `false`, allows continued tracing.
- `is_planar`: If `true`, applies planar surface projection `sqrt(max(0, k_hat · normal))`.
"""
mutable struct FourierBeamletPropagator{T <: Real, S <: AbstractShape{T}} <: AbstractFourierDetector{T}
    const shape::S
    hits::Union{
        Vector{AbstractDetectorHit},
        Vector{RayHit{T}},
        Vector{PolarizedRayHit{T}},
        Vector{GaussianBeamletHit{T}},
        Vector{AstigmaticGaussianBeamletHit{T}},
        Vector{<:AbstractDetectorHit}
    }
    stop::Bool
    is_planar::Bool
end

"""
    FourierBeamletPropagator(shape::AbstractShape{T}; stop::Bool = true, is_planar::Bool = false) where {T <: Real}
    FourierBeamletPropagator(shape::AbstractShape, stop::Bool)

Constructs a FourierBeamletPropagator with the specified geometric aperture shape.

Parameters:
- `shape`: Detector aperture surface shape of type `AbstractShape{T}`.
- `stop`: Flag indicating whether the detector absorbs/stops incident rays, defaults to `true`.
- `is_planar`: Flag indicating planar aperture projection geometry, defaults to `false`.

Returns:
- New instance of `FourierBeamletPropagator{T, typeof(shape)}`.
"""
FourierBeamletPropagator(shape::S; stop::Bool = true, is_planar::Bool = false) where {T <: Real, S <: AbstractShape{T}} =
    FourierBeamletPropagator{T, S}(shape, AbstractDetectorHit[], stop, is_planar)

FourierBeamletPropagator(shape::AbstractShape, stop::Bool) = FourierBeamletPropagator(shape; stop = stop)

"""
    FourierBeamletPropagator(edge_length::Real; stop::Bool = true, is_planar::Bool = false)
    FourierBeamletPropagator(edge_length::Real, stop::Bool)

Constructs a FourierBeamletPropagator with quadratic planar or cubic volume mesh geometry.

Parameters:
- `edge_length`: Edge length of detector aperture or volume in meters of type `Real`.
- `stop`: Ray absorption flag, defaults to `true`.
- `is_planar`: Flag indicating planar mesh geometry, defaults to `false`.

Returns:
- New `FourierBeamletPropagator` instance.
"""
function FourierBeamletPropagator(edge_length::Real; stop::Bool = true, is_planar::Bool = false)
    s = Float64(edge_length)
    if is_planar
        mesh = QuadraticFlatMesh(s); zrotate3d!(mesh, π)
        return FourierBeamletPropagator(mesh; stop = stop, is_planar = true)
    else
        mesh = CubeMesh(s); translate3d!(mesh, [-s / 2, -s / 2, -s / 2]); mesh.pos = zeros(Float64, 3)
        return FourierBeamletPropagator(mesh; stop = stop, is_planar = false)
    end
end
FourierBeamletPropagator(edge_length::Real, stop::Bool) = FourierBeamletPropagator(edge_length; stop = stop)

"""
    FourierBeamletPropagator(x_edge::Real, z_edge::Real; stop::Bool = true, is_planar::Bool = true)
    FourierBeamletPropagator(x_edge::Real, z_edge::Real, stop::Bool)

Constructs a FourierBeamletPropagator with rectangular planar mesh geometry.

Parameters:
- `x_edge`: Aperture dimension along x in meters of type `Real`.
- `z_edge`: Aperture dimension along z in meters of type `Real`.
- `stop`: Ray absorption flag, defaults to `true`.
- `is_planar`: Flag indicating planar mesh geometry, defaults to `true`.

Returns:
- New `FourierBeamletPropagator` instance.
"""
function FourierBeamletPropagator(x_edge::Real, z_edge::Real; stop::Bool = true, is_planar::Bool = true)
    mesh = RectangularFlatMesh(Float64(x_edge), Float64(z_edge)); zrotate3d!(mesh, π)
    return FourierBeamletPropagator(mesh; stop = stop, is_planar = is_planar)
end
FourierBeamletPropagator(x_edge::Real, z_edge::Real, stop::Bool) = FourierBeamletPropagator(x_edge, z_edge; stop = stop)

"""
    FourierBeamletPropagator(dx::Real, dy::Real, dz::Real; stop::Bool = true, is_planar::Bool = false)
    FourierBeamletPropagator(dx::Real, dy::Real, dz::Real, stop::Bool)

Constructs a FourierBeamletPropagator with cuboid volume mesh geometry.

Parameters:
- `dx`: Dimension along x in meters of type `Real`.
- `dy`: Dimension along y in meters of type `Real`.
- `dz`: Dimension along z in meters of type `Real`.
- `stop`: Ray absorption flag, defaults to `true`.
- `is_planar`: Planar geometry flag, defaults to `false`.

Returns:
- New `FourierBeamletPropagator` instance.
"""
function FourierBeamletPropagator(dx::Real, dy::Real, dz::Real; stop::Bool = true, is_planar::Bool = false)
    mesh = CuboidMesh(Float64(dx), Float64(dy), Float64(dz)); translate3d!(mesh, [-dx / 2, -dy / 2, -dz / 2]); mesh.pos = zeros(Float64, 3)
    return FourierBeamletPropagator(mesh; stop = stop, is_planar = is_planar)
end
FourierBeamletPropagator(dx::Real, dy::Real, dz::Real, stop::Bool) = FourierBeamletPropagator(dx, dy, dz; stop = stop)

"""
    FourierBeamletPropagator(object::AbstractObject; stop::Bool = true, is_planar::Bool = false)

Constructs a FourierBeamletPropagator extracting the geometric shape from an optical object.

Parameters:
- `object`: Optical component object of type `AbstractObject`.
- `stop`: Ray absorption flag, defaults to `true`.
- `is_planar`: Planar geometry flag, defaults to `false`.

Returns:
- New `FourierBeamletPropagator` instance.
"""
FourierBeamletPropagator(object::AbstractObject; stop::Bool = true, is_planar::Bool = false) =
    FourierBeamletPropagator(deepcopy(shape(object)); stop = stop, is_planar = is_planar)

# Detector interface
"""
    shape(d::FourierBeamletPropagator)

Returns the geometric aperture shape of the propagator detector.
"""
shape(d::FourierBeamletPropagator) = d.shape

"""
    hits(d::FourierBeamletPropagator)

Returns the container of captured detector hit records.
"""
hits(d::FourierBeamletPropagator) = d.hits

"""
    hits!(d::FourierBeamletPropagator, new)

Replaces the internal hits container of the propagator.
"""
hits!(d::FourierBeamletPropagator, new) = (d.hits = new)

"""
    stop(d::FourierBeamletPropagator)

Returns whether the detector absorbs incident optical rays.
"""
stop(d::FourierBeamletPropagator) = d.stop

"""
    Base.empty!(d::FourierBeamletPropagator)

Clears all stored hit records from the propagator detector.
"""
Base.empty!(d::FourierBeamletPropagator) = (empty!(d.hits); d)

"""
    Base.isempty(d::FourierBeamletPropagator)

Returns `true` if no hit records have been captured by the propagator.
"""
Base.isempty(d::FourierBeamletPropagator) = isempty(d.hits)

"""
    Base.length(d::FourierBeamletPropagator)

Returns the total count of captured detector hit records.
"""
Base.length(d::FourierBeamletPropagator) = length(d.hits)

const _ConcreteHit{T} = Union{RayHit{T}, PolarizedRayHit{T}, GaussianBeamletHit{T}, AstigmaticGaussianBeamletHit{T}}

const _PROPAGATOR_HIT_LOCK = Base.ReentrantLock()

@inline _push_or_widen!(h::Vector{H}, hit::H) where {H} = (push!(h, hit); h)

@inline function _push_or_widen!(h::Vector{T_el}, hit::H) where {T_el, H}
    if hit isa T_el
        push!(h, hit)
        return h
    else
        T_super = typejoin(T_el, H)
        new_hits = Vector{T_super}(undef, length(h) + 1)
        copyto!(new_hits, 1, h, 1, length(h))
        new_hits[end] = hit
        return new_hits
    end
end

@inline function _push_hit_unlocked!(d::FourierBeamletPropagator, hit::H) where {H}
    h = d.hits
    if isempty(h) && !(h isa Vector{H})
        d.hits = H[hit]
    else
        d.hits = _push_or_widen!(h, hit)
    end
    return d
end

"""
    Base.push!(d::FourierBeamletPropagator{T}, hit::AbstractDetectorHit) where {T}

Appends a detector hit record to the internal hits container of the propagator.
- Input: `d::FourierBeamletPropagator{T}`, `hit::AbstractDetectorHit`.
- Output: Mutated propagator instance `d`.
- Invariants: Preserves typed container when hit types match or dynamically widens type.
"""
function Base.push!(d::FourierBeamletPropagator{T}, hit::AbstractDetectorHit) where {T}
    lock(_PROPAGATOR_HIT_LOCK) do
        _push_hit_unlocked!(d, hit)
    end
    return d
end

"""
    Base.push!(d::FourierBeamletPropagator{T}, hit::H) where {T, H <: _ConcreteHit{T}}

Appends a concrete detector hit record to the internal hits container preserving zero heap allocation.
- Input: `d::FourierBeamletPropagator{T}`, `hit::H`.
- Output: Mutated propagator instance `d`.
- Invariants: Zero-allocation push when container has capacity (@allocated == 0).
"""
@inline function Base.push!(d::FourierBeamletPropagator{T}, hit::H) where {T, H <: _ConcreteHit{T}}
    h = d.hits
    if h isa Vector{H}
        push!(h, hit)
    elseif isempty(h)
        d.hits = H[hit]
    else
        invoke(Base.push!, Tuple{FourierBeamletPropagator{T}, AbstractDetectorHit}, d, hit)
    end
    return d
end

# ------------------------------------------------------------------------------
# interact3d methods
# ------------------------------------------------------------------------------

"""
    interact3d(::AbstractSystem, d::FourierBeamletPropagator, beam::Beam, ray::AbstractRay)

Intersects and records an optical ray hit on the propagator detector.
- Input: `sys::AbstractSystem`, `d::FourierBeamletPropagator`, `beam::Beam`, `ray::AbstractRay`.
- Output: `BeamInteraction` if `!d.stop`, or `nothing` if absorbing.
- Invariants: Appends hit to `d.hits`; returns unperturbed interaction when `!d.stop`.
"""
@inline _ray_hit(ray::PolarizedRay, opl) = PolarizedRayHit(ray, opl)
@inline _ray_hit(ray::AbstractRay, opl) = RayHit(ray, opl)

@inline function interact3d(::AbstractSystem, d::FourierBeamletPropagator, beam::Beam, ray::AbstractRay)
    opl = safe_finite(optical_path_length(beam))
    hit = _ray_hit(ray, opl)
    Base.push!(d, hit)
    return d.stop ? nothing : BeamInteraction(nothing, ray)
end

"""
    interact3d(::AbstractSystem, d::FourierBeamletPropagator{T}, g::GaussianBeamlet{R}, id::Int = 1) where {T <: Real, R}

Intersects and records a circular Gaussian beamlet hit on the propagator detector.
- Input: `sys::AbstractSystem`, `d::FourierBeamletPropagator{T}`, `g::GaussianBeamlet{R}`, `id::Int`.
- Output: `GaussianBeamletInteraction{R}` if `!d.stop`, or `nothing` if absorbing.
- Invariants: Appends GaussianBeamletHit to `d.hits`.
"""
function interact3d(::AbstractSystem, d::FourierBeamletPropagator{T}, g::GaussianBeamlet{R}, id::Int = 1) where {T <: Real, R}
    l_tot, l_cur = length(g), length(g.chief.rays[id])
    l0 = safe_finite(l_tot - l_cur, zero(R))
    ray = g.chief.rays[id]
    p0, d0 = position(ray), direction(ray)
    inter = intersection(ray)
    sqrt_proj = isnothing(inter) ? one(R) : sqrt(abs(dot(d0, normal3d(inter))))
    l_eval = safe_finite(length(g), zero(R))
    w_at_detector, _, _, _ = gauss_parameters(g, l_eval)
    w_max = safe_pos(w_at_detector, g.w0)
    Base.push!(d, GaussianBeamletHit(g, l0, id, p0, d0, sqrt_proj, w_max))
    return d.stop ? nothing : GaussianBeamletInteraction{R}(
        BeamInteraction(nothing, g.chief.rays[id]),
        BeamInteraction(nothing, g.waist.rays[id]),
        BeamInteraction(nothing, g.divergence.rays[id]),
    )
end

"""
    interact3d(::AbstractSystem, d::FourierBeamletPropagator{T}, agb::AstigmaticGaussianBeamlet{R}, id::Int = 1) where {T <: Real, R}

Intersects and records an astigmatic Gaussian beamlet hit on the propagator detector.
- Input: `sys::AbstractSystem`, `d::FourierBeamletPropagator{T}`, `agb::AstigmaticGaussianBeamlet{R}`, `id::Int`.
- Output: `AstigmaticGaussianBeamletInteraction{R}` if `!d.stop`, or `nothing` if absorbing.
- Invariants: Appends AstigmaticGaussianBeamletHit to `d.hits`.
"""
function interact3d(::AbstractSystem, d::FourierBeamletPropagator{T}, agb::AstigmaticGaussianBeamlet{R}, id::Int = 1) where {T <: Real, R}
    l_tot, l_cur = length(agb), length(agb.c.rays[id])
    l0 = safe_finite(l_tot - l_cur, zero(R))
    chief = rays(agb.c)[id]
    p0, d0, k0 = position(chief), direction(chief), 2π / wavelength(chief)
    inter = intersection(chief)
    sqrt_proj = isnothing(inter) ? one(R) : sqrt(abs(dot(d0, normal3d(inter))))
    h1, u1, h2, u2, _ = parabasal_ray_parameters(agb, p0, id)
    p0n, in_ = point_on_beam(agb, 0.0)
    dirn = direction(rays(agb.c)[in_])
    h1n, _, h2n, _, _ = parabasal_ray_parameters(agb, p0n, in_)
    area_ref = _pseudo_cross2d(h1n, h2n, dirn)
    E_vec = polarization(rays(agb.c)[in_])
    max_idx = argmax(abs.(E_vec))
    E_ref_amp = Complex{R}(norm(E_vec) * cis(angle(E_vec[max_idx])))
    p_parent = agb.parent
    l_parent = isnothing(p_parent) ? 0.0 : length(p_parent)
    opl_parent = isnothing(p_parent) ? 0.0 : optical_path_length(p_parent)
    Δl = opl_parent - l_parent
    for j in 1:(id - 1)
        ray_j = rays(agb.c)[j]
        Δl += optical_path_length(ray_j) - length(ray_j)
    end
    l_eval = safe_finite(length(chief), zero(length(chief)))
    H1, H2 = h1 + l_eval * u1, h2 + l_eval * u2
    w_max = safe_pos(max(norm(H1), norm(H2)), max(norm(h1), norm(h2)))
    n_eff = refractive_index(agb, id)
    Base.push!(d, AstigmaticGaussianBeamletHit(agb, l0, id, p0, d0, h1, u1, h2, u2,
        area_ref, k0, Δl, n_eff, E_ref_amp, sqrt_proj, w_max))
    return d.stop ? nothing : AstigmaticGaussianBeamletInteraction{R}(
        BeamInteraction(nothing, rays(agb.c)[id]),
        BeamInteraction(nothing, rays(agb.wxp)[id]),
        BeamInteraction(nothing, rays(agb.wxm)[id]),
        BeamInteraction(nothing, rays(agb.wyp)[id]),
        BeamInteraction(nothing, rays(agb.wym)[id]),
        BeamInteraction(nothing, rays(agb.dxp)[id]),
        BeamInteraction(nothing, rays(agb.dxm)[id]),
        BeamInteraction(nothing, rays(agb.dyp)[id]),
        BeamInteraction(nothing, rays(agb.dym)[id]),
    )
end

"""
    interact3d(d::FourierBeamletPropagator, ray::AbstractRay)
    interact3d(sys::AbstractSystem, d::FourierBeamletPropagator, ray::AbstractRay)
    interact3d(d::FourierBeamletPropagator, beam::Beam)
    interact3d(d::FourierBeamletPropagator, rays::AbstractVector{<:AbstractRay})
    interact3d(sys::AbstractSystem, d::FourierBeamletPropagator, beam::Beam)
    interact3d(sys::AbstractSystem, d::FourierBeamletPropagator, rays::AbstractVector{<:AbstractRay})
    interact3d(d::FourierBeamletPropagator, g::GaussianBeamlet, id::Int = 1)
    interact3d(d::FourierBeamletPropagator, agb::AstigmaticGaussianBeamlet, id::Int = 1)

Propagates optical rays, beams, and beamlets into the detector and registers interaction hits.

Parameters:
- `sys`: Optical system context of type `AbstractSystem`.
- `d`: Propagator detector instance of type `FourierBeamletPropagator`.
- `ray` / `beam` / `rays` / `g` / `agb`: Incident optical entity to interact with the detector.
- `id`: Beamlet ray index, defaults to 1.

Returns:
- `nothing` if detector stops rays, or the transmitted optical entity.
"""
interact3d(d::FourierBeamletPropagator, ray::AbstractRay) = interact3d(StaticSystem((d,)), d, Beam(ray), ray)
interact3d(sys::AbstractSystem, d::FourierBeamletPropagator, ray::AbstractRay) = interact3d(sys, d, Beam(ray), ray)
interact3d(d::FourierBeamletPropagator, beam::Beam) = (for r in rays(beam) interact3d(d, r) end; d.stop ? nothing : beam)
interact3d(d::FourierBeamletPropagator, rays::AbstractVector{<:AbstractRay}) = (for r in rays interact3d(d, r) end; d.stop ? nothing : rays)
interact3d(sys::AbstractSystem, d::FourierBeamletPropagator, beam::Beam) = interact3d(d, beam)
interact3d(sys::AbstractSystem, d::FourierBeamletPropagator, rays::AbstractVector{<:AbstractRay}) = interact3d(d, rays)
interact3d(d::FourierBeamletPropagator, g::GaussianBeamlet, id::Int = 1) = interact3d(StaticSystem((d,)), d, g, id)
interact3d(d::FourierBeamletPropagator, agb::AstigmaticGaussianBeamlet, id::Int = 1) = interact3d(StaticSystem((d,)), d, agb, id)

# ------------------------------------------------------------------------------
# Fourier Spectrum Generation & Hit Decomposition Subroutines
# ------------------------------------------------------------------------------

"""
    _evaluate_wave_parameters(hit::AbstractDetectorHit, ::Type{T}) where {T <: Real}

Evaluates optical wavelength, wavenumber magnitude, and angular frequency from a detector hit.

Parameters:
- `hit`: Detector hit record providing wavelength information.
- `T`: Floating-point numeric type for evaluated optical parameters.

Returns:
- A 3-tuple containing optical wavelength in meters, wavenumber magnitude in rad/m, and angular frequency in rad/s.

Contracts & Invariants:
- Evaluated parameters are strictly positive physical quantities.
"""
@inline function _evaluate_wave_parameters(hit::AbstractDetectorHit, ::Type{T}) where {T <: Real}
    wl = T(wavelength(hit)) # Centralized evaluation of wavelength in meters converted to precision T
    k0_mag = T(2π) / wl     # Fundamental wavenumber magnitude k0 = 2π / λ in rad/m
    omega = T(2π * c_vacuum) / wl # Optical angular frequency ω = 2π * c / λ in rad/s
    return (wl, k0_mag, omega)
end

"""
    _extract_hit_position(hit::AbstractDetectorHit, ::Type{T}) where {T <: Real}

Extracts the physical position vector of a detector hit with finite fallback.

Parameters:
- `hit`: Detector hit record providing hit coordinates or origin position.
- `T`: Floating-point numeric type for evaluated coordinate components.

Returns:
- A 3D static coordinate vector of type `SVector{3, T}` in meters.

Contracts & Invariants:
- All returned vector components are finite numbers.
"""
@inline function _extract_hit_position(hit::AbstractDetectorHit, ::Type{T}) where {T <: Real}
    hp = hit_point(hit)
    p_orig = position(hit)
    return SVector{3, T}(
        T(safe_finite(hp[1], p_orig[1])),
        T(safe_finite(hp[2], p_orig[2])),
        T(safe_finite(hp[3], p_orig[3]))
    )
end

"""
    _extract_hit_direction(hit::AbstractDetectorHit, ::Type{T}) where {T <: Real}

Extracts and normalizes the propagation direction vector from a detector hit.

Parameters:
- `hit`: Detector hit record providing propagation direction.
- `T`: Floating-point numeric type for evaluated vector components.

Returns:
- A normalized 3D static direction vector of type `SVector{3, T}`.

Contracts & Invariants:
- The returned vector has unit Euclidean norm.
"""
@inline function _extract_hit_direction(hit::AbstractDetectorHit, ::Type{T}) where {T <: Real}
    d_raw = direction(hit)
    return normalize(SVector{3, T}(T(d_raw[1]), T(d_raw[2]), T(d_raw[3]))) # Normalized propagation direction vector
end

"""
    _deduce_fourier_beamlet_sampling(sigma_k::T, L_box::T, dk::Union{Nothing, Real}) where {T}

Computes the spectral sampling step size and half-mode index count for discrete Fourier beamlet expansion.

Parameters:
- `sigma_k`: Spatial spectral bandwidth in rad/m of type `T`.
- `L_box`: Spatial box aperture dimension in meters of type `T`.
- `dk`: User-specified spectral step size in rad/m, or nothing for adaptive resolution.

Returns:
- A 2-tuple containing the spectral step size of type `T` and the integer half-mode index count.

Contracts & Invariants:
- The spectral step size is strictly positive.
- The mode count is a positive integer clamped within computational limits.
"""
function _deduce_fourier_beamlet_sampling(sigma_k::T, L_box::T, dk::Union{Nothing, Real}) where {T <: Real}
    if !isnothing(dk)
        dk_step = T(dk) # User-specified spectral step size converted to precision T
    else
        dk_box = T(2π) / L_box # Fundamental spatial-frequency box resolution in rad/m
        M_est = ceil(Int, 4 * sigma_k / dk_box) # Estimated half-mode count to cover 4*sigma_k bandwidth
        M_clamped = clamp(M_est, 1, 32) # Clamp mode half-count between 1 and 32 to bound computational complexity
        dk_step = (M_clamped == M_est) ? dk_box : (4 * sigma_k / T(M_clamped)) # Adjusted step size if mode count is clamped
    end
    M = max(1, ceil(Int, 4 * sigma_k / dk_step)) # Clamped half-mode sampling index bound
    return (dk_step, M)
end

function _deduce_fourier_beamlet_sampling(sigma_k::Real, L_box::Real, dk::Union{Nothing, Real})
    T = promote_type(typeof(float(sigma_k)), typeof(float(L_box)))
    return _deduce_fourier_beamlet_sampling(T(sigma_k), T(L_box), dk)
end


@inline function _extract_hit_ray_length(hit::GaussianBeamletHit, ::Type{T}) where {T <: Real}
    return (1 <= hit.id <= length(hit.gauss.chief.rays)) ? T(length(hit.gauss.chief.rays[hit.id])) : zero(T)
end
@inline function _extract_hit_ray_length(hit::AstigmaticGaussianBeamletHit, ::Type{T}) where {T <: Real}
    return (1 <= hit.id <= length(hit.agb.c.rays)) ? T(length(hit.agb.c.rays[hit.id])) : zero(T)
end
@inline _extract_hit_ray_length(hit, ::Type{T}) where {T <: Real} = zero(T)

"""
    _evaluate_hit_eval_length(hit, r_hit::SVector{3, T})::T where {T <: Real}

Determines the optical propagation evaluation distance for a beamlet hit.

Parameters:
- `hit`: Beamlet detector hit record.
- `r_hit`: Incident hit position 3-vector in meters of type `SVector{3, T}`.

Returns:
- Cumulative propagation path length in meters of type `T`.
"""
@inline function _evaluate_hit_eval_length(hit, r_hit::SVector{3, T})::T where {T <: Real}
    l_ray = _extract_hit_ray_length(hit, T)
    l_eval = safe_finite(T(hit.l0), zero(T)) + safe_finite(l_ray, zero(T))
    if iszero(l_eval)
        l_geom = norm(r_hit - position(hit))
        l_eval = safe_pos(T(l_geom), l_eval)
    end
    return l_eval
end

"""
    _beamlet_defocus_pair(hit::GaussianBeamletHit, k0_mag::T, r_hit::SVector{3, T}, w01::T, w02::T) where {T <: Real}
    _beamlet_defocus_pair(hit::AstigmaticGaussianBeamletHit, k0_mag::T, r_hit::SVector{3, T}, w01::T, w02::T) where {T <: Real}
    _beamlet_defocus_pair(hit, k0_mag::T, r_hit::SVector{3, T}, w01::T, w02::T) where {T <: Real}

Evaluates transverse defocus distances along primary axes for a beamlet hit.

# Arguments
- `hit`: Beamlet hit record.
- `k0_mag::T`: Wavenumber magnitude in rad/m of type `T`.
- `r_hit::SVector{3, T}`: Incident hit position vector in meters of type `SVector{3, T}`.
- `w01::T`: Primary waist radius in meters of type `T`.
- `w02::T`: Secondary waist radius in meters of type `T`.

# Returns
- Tuple `(dz1, dz2)` of defocus distances in meters of type `Tuple{T, T}`.
"""
@inline function _beamlet_defocus_pair(hit::GaussianBeamletHit, k0_mag::T, r_hit::SVector{3, T}, w01::T, w02::T) where {T <: Real}
    wl = T(2π) / k0_mag
    w, curv, _, _ = gauss_parameters(hit.gauss, _evaluate_hit_eval_length(hit, r_hit))
    dz = _beamlet_defocus_distance(safe_finite(T(curv), zero(T)), safe_pos(T(w), w01), wl)
    return (dz, dz)
end

@inline function _beamlet_defocus_pair(hit::AstigmaticGaussianBeamletHit, k0_mag::T, r_hit::SVector{3, T}, w01::T, w02::T) where {T <: Real}
    wl = T(2π) / k0_mag
    w1, w2, R1, R2, _, _, _ = gauss_parameters(hit.agb, _evaluate_hit_eval_length(hit, r_hit))
    return (_beamlet_defocus_distance(safe_inv(T(R1), zero(T)), safe_pos(T(w1), w01), wl),
            _beamlet_defocus_distance(safe_inv(T(R2), zero(T)), safe_pos(T(w2), w02), wl))
end

@inline _beamlet_defocus_pair(hit, k0_mag::T, r_hit::SVector{3, T}, w01::T, w02::T) where {T <: Real} = (zero(T), zero(T))

"""
    _sample_beamlet_modes!(pts, hit, dir, w01, w02, sigma_k, dk_step, M, t1, t2, E0, k0_mag, omega, r_hit, opl_val)

Samples transverse Gaussian beamlet Fourier modes on a two-dimensional grid and appends the resulting plane-wave modes to the output spectrum.

Parameters:
- `pts`: Output vector of Fourier k-points to be populated in-place.
- `hit`: Beamlet detector hit record.
- `dir`: Unit direction vector of the beamlet propagation axis.
- `w01`: First beam waist radius in meters along axis `t1`.
- `w02`: Second beam waist radius in meters along axis `t2`.
- `sigma_k`: Spatial spectral bandwidth in rad/m.
- `dk_step`: Transverse spectral sampling step in rad/m.
- `M`: Maximum mode sampling index bound.
- `t1`: First unit transverse polarization axis.
- `t2`: Second unit transverse polarization axis.
- `E0`: Complex reference electric field amplitude vector.
- `k0_mag`: Central wavenumber magnitude in rad/m.
- `omega`: Optical angular frequency in rad/s.
- `r_hit`: Position vector of the detector hit point in meters.
- `opl_val`: Accumulated optical path length in meters of type `T`.

Returns:
- Mutated output vector `pts` containing the appended Fourier k-points.

Contracts & Invariants:
- All generated k-points satisfy electromagnetic field transversality.
- Spectral weighting preserves Gaussian beamlet envelope normalization.
"""
function _sample_beamlet_modes!(
    pts::Vector{FourierKPoint{T}}, hit, dir::SVector{3, T}, w01::T, w02::T, sigma_k::T, dk_step::T,
    M::Int, t1::SVector{3, T}, t2::SVector{3, T}, E0::SVector{3, Complex{T}},
    k0_mag::T, omega::T, r_hit::SVector{3, T}, opl_val::T = zero(T)
) where {T <: Real}
    cutoff = 4 * sigma_k + T(1e-10)
    w01_sq_over_4 = w01^2 / T(4)
    w02_sq_over_4 = w02^2 / T(4)
    norm_factor = (sqrt(w01 * w02) * sqrt(w01 * w02) * dk_step^2) / (4 * T(π))
    dz1, dz2 = _beamlet_defocus_pair(hit, k0_mag, r_hit, w01, w02)
    inv_2k0 = inv(T(2) * k0_mag)
    for m1 in -M:M
        dk1 = m1 * dk_step
        for m2 in -M:M
            dk2 = m2 * dk_step
            dk_perp_sq = dk1^2 + dk2^2
            if sqrt(dk_perp_sq) <= cutoff
                k_vec = sqrt(max(zero(T), k0_mag^2 - dk_perp_sq)) * dir + dk1 * t1 + dk2 * t2
                k_hat = normalize(k_vec)
                phase_defocus = cis(-(dz1 * dk1^2 + dz2 * dk2^2) * inv_2k0)
                phase = cis(k0_mag * opl_val - dot(k_vec, r_hit))
                g = exp(-(w01_sq_over_4 * dk1^2 + w02_sq_over_4 * dk2^2))
                E_trans = E0 - dot(k_hat, E0) * k_hat
                E_vec = (norm_factor * g * phase * phase_defocus) * E_trans
                push!(pts, FourierKPoint{T}(k_vec, omega, E_vec, T(g)))
            end
        end
    end
    return pts
end

@inline _sample_beamlet_modes!(pts::Vector{FourierKPoint{T}}, hit, dir::SVector{3, T}, w0::T, sigma_k::T, dk_step::T, M::Int, t1::SVector{3, T}, t2::SVector{3, T}, E0::SVector{3, Complex{T}}, k0_mag::T, omega::T, r_hit::SVector{3, T}, opl_val::T = zero(T)) where {T <: Real} =
    _sample_beamlet_modes!(pts, hit, dir, w0, w0, sigma_k, dk_step, M, t1, t2, E0, k0_mag, omega, r_hit, opl_val)

"""
    _apply_planar_obliquity!(pts::Vector{FourierKPoint{T}}, n_hat::SVector{3, T}) where {T <: Real}

Applies the Lambertian planar obliquity cosine factor `sqrt(max(0, |k_hat · n_hat|))` to plane-wave modes in-place.
"""
@inline function _apply_planar_obliquity!(pts::Vector{FourierKPoint{T}}, n_hat::SVector{3, T}) where {T <: Real}
    @inbounds for i in eachindex(pts)
        p = pts[i]
        k_hat = direction(p)
        cos_fac = sqrt(max(zero(T), dot(k_hat, n_hat)))
        pts[i] = FourierKPoint{T}(p.k4, cos_fac * p.E, p.weight)
    end
    return pts
end

@inline function _detector_normal_local(d::FourierBeamletPropagator{S, Sh}, n_hat::SVector{3, T}) where {S <: Real, Sh <: AbstractShape, T <: Real}
    if d.is_planar && d.shape isa AbstractShape
        o = orientation(d.shape)
        R_trans = SMatrix{3, 3, T, 9}(
            T(o[1, 1]), T(o[1, 2]), T(o[1, 3]),
            T(o[2, 1]), T(o[2, 2]), T(o[2, 3]),
            T(o[3, 1]), T(o[3, 2]), T(o[3, 3])
        )'
        return normalize(R_trans * n_hat)
    else
        return n_hat
    end
end
@inline _detector_normal_local(::Bool, n_hat::SVector{3, T}) where {T <: Real} = n_hat

"""
    _transform_to_detector_frame(d::FourierBeamletPropagator{S, Sh}, dir::SVector{3, T}, E0::SVector{3, Complex{T}}, r_hit::SVector{3, T}) where {S <: Real, Sh <: AbstractShape, T <: Real}

Transforms directional and spatial ray coordinates from the global reference frame into the detector local frame.

Parameters:
- `d`: Propagator instance of type `FourierBeamletPropagator{S, Sh}`.
- `dir`: Propagation unit direction 3-vector of type `SVector{3, T}`.
- `E0`: Complex electric field polarization 3-vector of type `SVector{3, Complex{T}}`.
- `r_hit`: Incident hit position 3-vector in meters of type `SVector{3, T}`.

Returns:
- Tuple `(dir_local, E0_local, r_local, s_poca)` containing transformed unit direction, electric field, relative hit position, and optical distance to POCA.

Contracts & Invariants:
- Preserves the norm of the propagation direction vector.
- Preserves transversality and magnitude of the electric field vector.
- Inlined and allocation-free (@allocated == 0).
"""
@inline function _transform_to_detector_frame(
    d::FourierBeamletPropagator{S, Sh},
    dir::SVector{3, T},
    E0::SVector{3, Complex{T}},
    r_hit::SVector{3, T}
) where {S <: Real, Sh <: AbstractShape, T <: Real}
    c = SVector{3, T}(position(d.shape))
    if d.is_planar
        o = orientation(d.shape)
        R_trans = SMatrix{3, 3, T, 9}(
            T(o[1, 1]), T(o[1, 2]), T(o[1, 3]),
            T(o[2, 1]), T(o[2, 2]), T(o[2, 3]),
            T(o[3, 1]), T(o[3, 2]), T(o[3, 3])
        )'
        r_rel = r_hit - c
        return (normalize(R_trans * dir), R_trans * E0, R_trans * r_rel, zero(T))
    else
        s_poca = dot(c - r_hit, dir)
        r_poca = r_hit + s_poca * dir
        r_rel = r_poca - c
        return (dir, E0, r_rel, s_poca)
    end
end

@inline _transform_to_detector_frame(::Bool, dir::SVector{3, T}, E0::SVector{3, Complex{T}}, r_hit::SVector{3, T}) where {T <: Real} = (dir, E0, r_hit, zero(T))


"""
    _hit_power(hit::PolarizedRayHit, ::Type{T}) where {T <: Real}
    _hit_power(hit::AbstractRayHit, ::Type{T}) where {T <: Real}
    _hit_power(hit, ::Type{T}) where {T <: Real}

Extracts the optical power of a ray hit in precision `T`.

# Arguments
- `hit`: Ray hit record.
- `::Type{T}`: Target numeric type.

# Returns
- Optical power of type `T`.
"""
@inline _hit_power(hit::PolarizedRayHit, ::Type{T}) where {T <: Real} = T(sum(abs2, hit.ray.E0))
@inline _hit_power(hit::AbstractRayHit, ::Type{T}) where {T <: Real} = hasproperty(hit, :power) ? T(hit.power) : one(T)
@inline _hit_power(hit, ::Type{T}) where {T <: Real} = hasproperty(hit, :power) ? T(hit.power) : one(T)

"""
    _hit_polarization_raw(hit::PolarizedRayHit, dir)
    _hit_polarization_raw(hit, dir)

Extracts the electric field polarization vector for a ray hit.

# Arguments
- `hit`: Ray hit record.
- `dir`: Propagation direction vector.

# Returns
- Polarization vector.
"""
@inline _hit_polarization_raw(hit::PolarizedRayHit, dir) = hit.ray.E0
@inline _hit_polarization_raw(hit, dir) = transverse_polarization(dir)

"""
    _hit_to_kpoints!(pts::Vector{FourierKPoint{T}}, hit::AbstractRayHit, L_box::T, n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing) where {T <: Real}

Converts an incident ray hit into a discrete delta Fourier k-point mode and appends it to the spectrum.

Parameters:
- `pts`: Output vector of Fourier k-points to be populated in-place.
- `hit`: Ray detector hit record of type `AbstractRayHit`.
- `L_box`: Spatial box aperture dimension in meters of type `T`.
- `n_hat`: Unit surface normal vector of the detector aperture.
- `d`: Propagator instance or planar flag of type `Union{FourierBeamletPropagator, Bool}`.
- `dk`: User-specified spectral step size in rad/m, or nothing.

Returns:
- Mutated output vector `pts` containing the appended Fourier k-point.

Contracts & Invariants:
- Preserves optical power, phase, and polarization state.
"""
function _hit_to_kpoints!(
    pts::Vector{FourierKPoint{T}}, hit::AbstractRayHit, L_box::T,
    n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing
) where {T <: Real}
    dir = _extract_hit_direction(hit, T)
    wl, k0_mag, omega = _evaluate_wave_parameters(hit, T)
    r_hit = _extract_hit_position(hit, T)
    pol_raw = _hit_polarization_raw(hit, dir)
    p_vec = SVector{3, Complex{T}}(Complex{T}(pol_raw[1]), Complex{T}(pol_raw[2]), Complex{T}(pol_raw[3]))
    p_norm = norm(p_vec)
    e_hat = p_norm > 0 ? (p_vec / p_norm) : SVector{3, Complex{T}}(transverse_polarization(dir))
    dir_loc, e_loc, r_loc, s_poca = _transform_to_detector_frame(d, dir, e_hat, r_hit)
    k0 = k0_mag * dir_loc
    opl_val = safe_finite(T(optical_path_length(hit)), zero(T)) + s_poca
    phase = cis(k0_mag * opl_val - dot(k0, r_loc))
    P_ray = max(zero(T), _hit_power(hit, T))
    amp = sqrt(P_ray) * L_box
    E = amp * phase * e_loc
    push!(pts, FourierKPoint{T}(k0, omega, E, one(T)))
    return pts
end

"""
    _beamlet_waist_and_axes(hit::GaussianBeamletHit, dir::SVector{3, T}, r_hit::SVector{3, T}, ::Type{T}) where {T <: Real}
    _beamlet_waist_and_axes(hit::AstigmaticGaussianBeamletHit, dir::SVector{3, T}, r_hit::SVector{3, T}, ::Type{T}) where {T <: Real}
    _beamlet_waist_and_axes(hit::AbstractBeamletHit, dir::SVector{3, T}, r_hit::SVector{3, T}, ::Type{T}) where {T <: Real}

Determines transverse waist radii and orthogonal transverse axis vectors for a beamlet hit.

# Arguments
- `hit`: Beamlet hit record.
- `dir`: Propagation direction unit vector of type `SVector{3, T}`.
- `r_hit`: Incident hit position vector in meters of type `SVector{3, T}`.
- `::Type{T}`: Target numeric type.

# Returns
- Tuple `(w01, w02, t1_v, t2_v)` of waist radii and orthonormal transverse vectors.
"""
@inline function _beamlet_waist_and_axes(hit::GaussianBeamletHit, dir::SVector{3, T}, r_hit::SVector{3, T}, ::Type{T}) where {T <: Real}
    w = safe_pos(T(hit.gauss.w0), T(hit.w_max))
    t1 = SVector{3, T}(transverse_polarization(dir))
    return (w, w, t1, normalize(cross(dir, t1)))
end

@inline function _beamlet_waist_and_axes(hit::AstigmaticGaussianBeamletHit, dir::SVector{3, T}, r_hit::SVector{3, T}, ::Type{T}) where {T <: Real}
    _, _, _, _, _, w01, w02 = gauss_parameters(hit.agb, _evaluate_hit_eval_length(hit, r_hit))
    x_axis = SVector{3, T}(one(T), zero(T), zero(T))
    proj = x_axis - dot(dir, x_axis) * dir
    t1 = norm(proj) > T(1e-6) ? normalize(proj) : SVector{3, T}(zero(T), zero(T), one(T))
    return (safe_pos(T(w01), T(hit.w_max)), safe_pos(T(w02), T(hit.w_max)), t1, normalize(cross(dir, t1)))
end

@inline function _beamlet_waist_and_axes(hit::AbstractBeamletHit, dir::SVector{3, T}, r_hit::SVector{3, T}, ::Type{T}) where {T <: Real}
    w = safe_pos(T(hit.w_max), T(2π / hit.k0))
    t1 = SVector{3, T}(transverse_polarization(dir))
    return (w, w, t1, normalize(cross(dir, t1)))
end

"""
    _hit_to_kpoints!(pts::Vector{FourierKPoint{T}}, hit::AbstractBeamletHit, L_box::T, n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing) where {T <: Real}

Decomposes an incident Gaussian beamlet hit into discrete Fourier k-points and appends them to the spectrum.

Parameters:
- `pts`: Output vector of Fourier k-points to be populated in-place.
- `hit`: Beamlet detector hit record of type `AbstractBeamletHit`.
- `L_box`: Spatial box aperture dimension in meters of type `T`.
- `n_hat`: Unit surface normal vector of the detector aperture.
- `d`: Propagator instance or planar flag of type `Union{FourierBeamletPropagator, Bool}`.
- `dk`: User-specified spectral step size in rad/m, or nothing for adaptive resolution.

Returns:
- Mutated output vector `pts` containing the appended Fourier k-points.

Contracts & Invariants:
- All generated k-points satisfy electromagnetic field transversality.
- Evaluates Gaussian transverse spectral distribution centered around the beamlet axis.
"""
function _hit_to_kpoints!(
    pts::Vector{FourierKPoint{T}}, hit::AbstractBeamletHit, L_box::T,
    n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing
) where {T <: Real}
    wl, k0_mag, omega = _evaluate_wave_parameters(hit, T)
    dir = _extract_hit_direction(hit, T)
    r_hit = _extract_hit_position(hit, T)

    pol = polarization(hit)
    E0_raw = SVector{3, Complex{T}}(Complex{T}(pol[1]), Complex{T}(pol[2]), Complex{T}(pol[3]))
    E0_init = norm(E0_raw) > 0 ? E0_raw : SVector{3, Complex{T}}(transverse_polarization(dir))

    dir_loc, E0, r_loc, s_poca = _transform_to_detector_frame(d, dir, E0_init, r_hit)

    w01, w02, t1, t2 = _beamlet_waist_and_axes(hit, dir_loc, r_hit, T)
    w01 = max(w01, T(1e-12))
    w02 = max(w02, T(1e-12))
    sigma_k = sqrt(T(2)) / min(w01, w02)

    dk_step, M = _deduce_fourier_beamlet_sampling(sigma_k, L_box, dk)

    opl_val = safe_finite(T(optical_path_length(hit)), zero(T)) + s_poca
    _sample_beamlet_modes!(pts, hit, dir_loc, w01, w02, sigma_k, dk_step, M, t1, t2, E0, k0_mag, omega, r_loc, opl_val)

    return pts
end

"""
    _hit_to_kpoints!(pts::Vector{FourierKPoint{T}}, hit::FourierKPoint, L_box::T, n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing) where {T <: Real}

Passes an existing Fourier k-point through to the spectrum with target type conversion.

Parameters:
- `pts`: Output vector of Fourier k-points to be populated in-place.
- `hit`: Existing Fourier k-point instance of type `FourierKPoint`.
- `L_box`: Spatial box aperture dimension in meters of type `T`.
- `n_hat`: Unit surface normal vector of the detector aperture.
- `d`: Propagator instance or planar flag of type `Union{FourierBeamletPropagator, Bool}`.
- `dk`: User-specified spectral step size in rad/m, or nothing.

Returns:
- Mutated output vector `pts` containing the appended Fourier k-point.

Contracts & Invariants:
- Preserves wavevector, frequency, electric field amplitude, and spectral weight.
- Converts all fields to precision `T`.
"""
function _hit_to_kpoints!(
    pts::Vector{FourierKPoint{T}}, hit::FourierKPoint, L_box::T,
    n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing
) where {T <: Real}
    k4_T = SVector{4, T}(T(hit.k4[1]), T(hit.k4[2]), T(hit.k4[3]), T(hit.k4[4])) # Four-wavevector converted to precision T
    E_T = SVector{3, Complex{T}}(Complex{T}(hit.E[1]), Complex{T}(hit.E[2]), Complex{T}(hit.E[3])) # Electric field converted to precision T
    push!(pts, FourierKPoint{T}(k4_T, E_T, T(hit.weight))) # Pass-through existing k-point
    return pts
end

"""
    _hit_to_kpoints!(pts::Vector{FourierKPoint{T}}, hit::AbstractDetectorHit, L_box::T, n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing) where {T <: Real}

Fallback handler for unsupported detector hit types.

Parameters:
- `pts`: Output vector of Fourier k-points.
- `hit`: Arbitrary detector hit record.
- `L_box`: Spatial box aperture dimension.
- `n_hat`: Unit surface normal vector.
- `d`: Propagator instance or planar flag of type `Union{FourierBeamletPropagator, Bool}`.
- `dk`: Spectral sampling step size parameter.

Returns:
- Mutated output vector `pts` unchanged.

Contracts & Invariants:
- Unrecognized hit types are safely ignored without raising runtime exceptions.
"""
function _hit_to_kpoints!(
    pts::Vector{FourierKPoint{T}}, hit::AbstractDetectorHit, L_box::T,
    n_hat::SVector{3, T}, d::Union{FourierBeamletPropagator, Bool}, dk::Union{Nothing, Real} = nothing
) where {T <: Real}
    return pts
end

"""
    generate_fourier_spectrum(d::FourierBeamletPropagator{T}, hits = d.hits;
                              L_box = 1.0, normal = SVector{3, T}(0, 1, 0), dk = nothing) where {T <: Real}

Generates the discrete Fourier angular spectrum `Vector{FourierKPoint{T}}` from captured hits.

Parameters:
- `d`: Propagator instance of type `FourierBeamletPropagator{T}`.
- `hits`: Captured detector hits to transform into Fourier k-points, defaulting to `d.hits`.
- `L_box`: Spatial box aperture dimension in meters of type `Real`, defaults to 1.0.
- `normal`: Surface normal vector of the detector aperture of type `AbstractVector`, defaults to `[0, 1, 0]`.
- `dk`: User-specified spectral sampling step size in rad/m of type `Union{Nothing, Real}`. If `nothing`, the step size is adaptively deduced from `L_box` and the beamlet spatial spectral bandwidth `sigma_k`, bounded to a maximum of 32 half-modes per transverse dimension. If a numeric value is provided, it explicitly defines the uniform transverse spectral sampling step `dk_step`.

Returns:
- Vector of discrete Fourier k-points of type `Vector{FourierKPoint{T}}`.

Contracts & Invariants:
- For simple rays, generates delta k-points along the ray propagation direction.
- For Gaussian beamlets, evaluates transverse Gaussian spectral modes enforcing transversality `k · E = 0`.
- If `d.is_planar` is true, applies obliquity projection scaling using the aperture normal.
- Preserves energy and phase anchoring relative to the hit position.
"""
function generate_fourier_spectrum(
    d::FourierBeamletPropagator{T},
    hits::AbstractVector = d.hits;
    L_box::Real = 1.0,
    normal::AbstractVector = SVector{3, T}(0, 1, 0),
    dk::Union{Nothing, Real} = nothing
) where {T <: Real}
    pts = FourierKPoint{T}[]
    isempty(hits) && return pts

    L_box_T = T(L_box) # Box size converted to precision T
    n_hat = normalize(SVector{3, T}(T(normal[1]), T(normal[2]), T(normal[3]))) # Unit aperture normal vector

    for hit in hits
        _hit_to_kpoints!(pts, hit, L_box_T, n_hat, d, dk)
    end

    if d.is_planar
        _apply_planar_obliquity!(pts, _detector_normal_local(d, n_hat))
    end

    return pts
end

# ==============================================================================
# Layer 3: Fourier Gridding & Field Synthesis (synthesize_field)
# ==============================================================================

using NonuniformFFTs

"""
    AbstractGrid{T <: Real}

Abstract supertype for spatial detector grids.
"""
abstract type AbstractGrid{T <: Real} end

"""
    SpatialGrid{T <: Real, D} <: AbstractGrid{T}

Spatial grid defined over D coordinate ranges.
Fields:
- `ranges`: Coordinate LinRanges for each dimension.
"""
struct SpatialGrid{T <: Real, D} <: AbstractGrid{T}
    ranges::NTuple{D, LinRange{T, Int}}

    @inline SpatialGrid{T, D}(ranges::NTuple{D, LinRange{T, Int}}) where {T <: Real, D} = new{T, D}(ranges)
end

@inline function SpatialGrid(ranges::NTuple{D, AbstractRange}) where {D}
    T_range = promote_type(map(eltype, ranges)...)
    T = T_range <: Real ? T_range : Float64
    lin_ranges = map(r -> LinRange{T, Int}(T(first(r)), T(last(r)), length(r)), ranges)
    return SpatialGrid{T, D}(lin_ranges)
end

@inline SpatialGrid(ranges::AbstractRange...) = SpatialGrid(ranges)

@inline function SpatialGrid{T, D}(ranges::NTuple{D, AbstractRange}) where {T <: Real, D}
    lin_ranges = map(r -> LinRange{T, Int}(T(first(r)), T(last(r)), length(r)), ranges)
    return SpatialGrid{T, D}(lin_ranges)
end

@inline SpatialGrid{T}(ranges::NTuple{D, AbstractRange}) where {T <: Real, D} = SpatialGrid{T, D}(ranges)

"""
    _grid_step_and_center(r::LinRange{T}) where {T <: Real}

Computes coordinate spacing step and central coordinate of a 1D spatial grid range.

Parameters:
- `r`: Linear coordinate range of type `LinRange{T}`.

Returns:
- 2-tuple containing step size and center coordinate of type `(T, T)`.
"""
@inline _grid_step_and_center(r::LinRange{T}) where {T <: Real} = (T(step(r)), first(r) + (length(r) ÷ 2) * T(step(r)))

"""
    _synthesize_field_nufft!(E_out::AbstractArray{SVector{3, Complex{T}}, 1}, pts::AbstractVector{<:FourierKPoint{T}}, grid::SpatialGrid{T, 1}) where {T}

In-place synthesis of the 1D complex vector electric field on the spatial detector grid.

Parameters:
- `E_out`: Destination array for complex 3D electric field vectors in V/m.
- `pts`: Discrete Fourier k-point spectrum of type AbstractVector{<:FourierKPoint{T}}.
- `grid`: One-dimensional spatial detector grid of type SpatialGrid{T, 1}.

Output:
- Mutated destination array `E_out` containing the synthesized electric field in V/m.

Contracts & Invariants:
- Grid dimension length must be at least 4.
- In-place array dimension matches the grid dimension.
"""
function _synthesize_field_nufft!(E_out::AbstractArray{SVector{3, Complex{T}}, 1}, pts::AbstractVector{<:FourierKPoint{T}}, grid::SpatialGrid{T, 1}, ::Union{FourierBeamletPropagator, Nothing} = nothing) where {T}
    r1 = grid.ranges[1]
    N1 = length(r1)
    if N1 < 4
        throw(ArgumentError("SpatialGrid dimensions must be >= 4 for NUFFT field synthesis"))
    end
    if isempty(pts)
        fill!(E_out, zero(SVector{3, Complex{T}}))
        return E_out
    end
    dx1, xc1 = _grid_step_and_center(r1)
    n_pts = length(pts)
    pts1 = Vector{T}(undef, n_pts)
    v1 = Vector{Complex{T}}(undef, n_pts)
    v2 = Vector{Complex{T}}(undef, n_pts)
    v3 = Vector{Complex{T}}(undef, n_pts)
    @inbounds for j in 1:n_pts
        p = pts[j]
        pts1[j] = mod(-p.k[1] * dx1, T(2π))
        cp = cis(p.k[1] * xc1)
        v1[j] = p.E[1] * cp
        v2[j] = p.E[2] * cp
        v3[j] = p.E[3] * cp
    end
    plan = PlanNUFFT(Complex{T}, (N1,); ntransforms = Val(3), fftshift = true)
    set_points!(plan, pts1)
    u1 = Array{Complex{T}, 1}(undef, N1)
    u2 = Array{Complex{T}, 1}(undef, N1)
    u3 = Array{Complex{T}, 1}(undef, N1)
    exec_type1!((u1, u2, u3), plan, (v1, v2, v3))
    @inbounds for i1 in 1:N1
        E_out[i1] = SVector{3, Complex{T}}(u1[i1], u2[i1], u3[i1])
    end
    return E_out
end

"""
    _detector_local_axes(d::FourierBeamletPropagator{S, Sh}, ::Type{T} = S) where {S <: Real, Sh <: AbstractShape, T <: Real}

Extracts the orthonormal transverse basis vectors `(ex, ez)` for the detector aperture or canonical projection axes.

Parameters:
- `d`: The Fourier beamlet propagator detector instance.
- `::Type{T}`: Target numeric type for coordinate axes, defaulting to `S`.

Output:
- A 2-tuple `(ex, ez)` of orthonormal unit vectors of type `NTuple{2, SVector{3, T}}` defining transverse spatial coordinates.

Contracts & Invariants:
- `ex` and `ez` are orthonormal unit vectors spanning the transverse detector or grid plane.
"""
@inline function _detector_local_axes(d::FourierBeamletPropagator{S, Sh}, ::Type{T} = S) where {S <: Real, Sh <: AbstractShape, T <: Real}
    if d.is_planar && d.shape isa AbstractShape
        o = orientation(d.shape)
        ex = SVector{3, T}(-T(o[1, 1]), -T(o[2, 1]), -T(o[3, 1]))
        ez = SVector{3, T}(T(o[1, 3]), T(o[2, 3]), T(o[3, 3]))
        return (ex, ez)
    else
        return (SVector{3, T}(1, 0, 0), SVector{3, T}(0, 1, 0))
    end
end

"""
    _synthesize_field_nufft!(E_out::AbstractArray{SVector{3, Complex{T}}, 2}, pts::AbstractVector{<:FourierKPoint{T}}, grid::SpatialGrid{T, 2}, d::Union{FourierBeamletPropagator, Nothing} = nothing) where {T}

In-place synthesis of the 2D complex vector electric field on the spatial detector grid.

Parameters:
- `E_out`: Destination array for complex 3D electric field vectors in V/m.
- `pts`: Discrete Fourier k-point spectrum of type AbstractVector{<:FourierKPoint{T}}.
- `grid`: Two-dimensional spatial detector grid of type SpatialGrid{T, 2}.
- `d`: Optional detector instance for aperture orientation alignment.

Output:
- Mutated destination array `E_out` containing the synthesized electric field in V/m.

Contracts & Invariants:
- Grid dimension lengths must be at least 4.
- In-place array dimensions match the grid dimensions.
"""
function _synthesize_field_nufft!(
    E_out::AbstractArray{SVector{3, Complex{T}}, 2},
    pts::AbstractVector{<:FourierKPoint{T}},
    grid::SpatialGrid{T, 2},
    d::Union{FourierBeamletPropagator, Nothing} = nothing
) where {T}
    r1, r2 = grid.ranges[1], grid.ranges[2]
    N1, N2 = length(r1), length(r2)
    if N1 < 4 || N2 < 4
        throw(ArgumentError("SpatialGrid dimensions must be >= 4 for NUFFT field synthesis"))
    end
    if isempty(pts)
        fill!(E_out, zero(SVector{3, Complex{T}}))
        return E_out
    end
    dx1, xc1 = _grid_step_and_center(r1)
    dx2, xc2 = _grid_step_and_center(r2)
    ex, ez = isnothing(d) ? (SVector{3, T}(1, 0, 0), SVector{3, T}(0, 1, 0)) : _detector_local_axes(d, T)
    n_pts = length(pts)
    pts1 = Vector{T}(undef, n_pts)
    pts2 = Vector{T}(undef, n_pts)
    v1 = Vector{Complex{T}}(undef, n_pts)
    v2 = Vector{Complex{T}}(undef, n_pts)
    v3 = Vector{Complex{T}}(undef, n_pts)
    @inbounds for j in 1:n_pts
        p = pts[j]
        ku = dot(p.k, ex)
        kv = dot(p.k, ez)
        pts1[j] = mod(-ku * dx1, T(2π))
        pts2[j] = mod(-kv * dx2, T(2π))
        cp = cis(ku * xc1 + kv * xc2)
        v1[j] = p.E[1] * cp
        v2[j] = p.E[2] * cp
        v3[j] = p.E[3] * cp
    end
    plan = PlanNUFFT(Complex{T}, (N1, N2); ntransforms = Val(3), fftshift = true)
    set_points!(plan, (pts1, pts2))
    u1 = Array{Complex{T}, 2}(undef, N1, N2)
    u2 = Array{Complex{T}, 2}(undef, N1, N2)
    u3 = Array{Complex{T}, 2}(undef, N1, N2)
    exec_type1!((u1, u2, u3), plan, (v1, v2, v3))
    @inbounds for i2 in 1:N2
        for i1 in 1:N1
            E_out[i1, i2] = SVector{3, Complex{T}}(u1[i1, i2], u2[i1, i2], u3[i1, i2])
        end
    end
    return E_out
end

"""
    _synthesize_field_nufft!(E_out::AbstractArray{SVector{3, Complex{T}}, 3}, pts::AbstractVector{<:FourierKPoint{T}}, grid::SpatialGrid{T, 3}) where {T}

In-place synthesis of the 3D complex vector electric field on the spatial detector grid.

Parameters:
- `E_out`: Destination array for complex 3D electric field vectors in V/m.
- `pts`: Discrete Fourier k-point spectrum of type AbstractVector{<:FourierKPoint{T}}.
- `grid`: Three-dimensional spatial detector grid of type SpatialGrid{T, 3}.

Output:
- Mutated destination array `E_out` containing the synthesized electric field in V/m.

Contracts & Invariants:
- Grid dimension lengths must be at least 4.
- In-place array dimensions match the grid dimensions.
"""
function _synthesize_field_nufft!(E_out::AbstractArray{SVector{3, Complex{T}}, 3}, pts::AbstractVector{<:FourierKPoint{T}}, grid::SpatialGrid{T, 3}, ::Union{FourierBeamletPropagator, Nothing} = nothing) where {T}
    r1, r2, r3 = grid.ranges[1], grid.ranges[2], grid.ranges[3]
    N1, N2, N3 = length(r1), length(r2), length(r3)
    if N1 < 4 || N2 < 4 || N3 < 4
        throw(ArgumentError("SpatialGrid dimensions must be >= 4 for NUFFT field synthesis"))
    end
    if isempty(pts)
        fill!(E_out, zero(SVector{3, Complex{T}}))
        return E_out
    end
    dx1, xc1 = _grid_step_and_center(r1)
    dx2, xc2 = _grid_step_and_center(r2)
    dx3, xc3 = _grid_step_and_center(r3)
    n_pts = length(pts)
    pts1 = Vector{T}(undef, n_pts)
    pts2 = Vector{T}(undef, n_pts)
    pts3 = Vector{T}(undef, n_pts)
    v1 = Vector{Complex{T}}(undef, n_pts)
    v2 = Vector{Complex{T}}(undef, n_pts)
    v3 = Vector{Complex{T}}(undef, n_pts)
    @inbounds for j in 1:n_pts
        p = pts[j]
        pts1[j] = mod(-p.k[1] * dx1, T(2π))
        pts2[j] = mod(-p.k[2] * dx2, T(2π))
        pts3[j] = mod(-p.k[3] * dx3, T(2π))
        cp = cis(p.k[1] * xc1 + p.k[2] * xc2 + p.k[3] * xc3)
        v1[j] = p.E[1] * cp
        v2[j] = p.E[2] * cp
        v3[j] = p.E[3] * cp
    end
    plan = PlanNUFFT(Complex{T}, (N1, N2, N3); ntransforms = Val(3), fftshift = true)
    set_points!(plan, (pts1, pts2, pts3))
    u1 = Array{Complex{T}, 3}(undef, N1, N2, N3)
    u2 = Array{Complex{T}, 3}(undef, N1, N2, N3)
    u3 = Array{Complex{T}, 3}(undef, N1, N2, N3)
    exec_type1!((u1, u2, u3), plan, (v1, v2, v3))
    @inbounds for i3 in 1:N3
        for i2 in 1:N2
            for i1 in 1:N1
                E_out[i1, i2, i3] = SVector{3, Complex{T}}(u1[i1, i2, i3], u2[i1, i2, i3], u3[i1, i2, i3])
            end
        end
    end
    return E_out
end

"""
    _extract_or_generate_kpoints(d::FourierBeamletPropagator{T}) where {T}

Centralizes the extraction or generation of Fourier k-points from propagator hits.

Parameters:
- `d`: Fourier beamlet propagator instance containing captured hits.

Returns:
- Vector of Fourier k-points of type `Vector{FourierKPoint{T}}`.

Contracts & Invariants:
- Output vector elements are strictly type-stable with element type `FourierKPoint{T}`.
- If hits already consist of `FourierKPoint`, avoids redundant regeneration.
"""
function _extract_or_generate_kpoints(d::FourierBeamletPropagator{T})::Vector{FourierKPoint{T}} where {T}
    if d.hits isa Vector{FourierKPoint{T}}
        return d.hits
    elseif !isempty(d.hits) && all(h -> h isa FourierKPoint, d.hits)
        return FourierKPoint{T}[FourierKPoint{T}(p.k4, p.E, p.weight) for p in d.hits]
    else
        return generate_fourier_spectrum(d)
    end
end

"""
    _extract_or_generate_kpoints(d::FourierBeamletPropagator, ::Type{T}) where {T <: Real}

Centralizes the extraction or generation of Fourier k-points converted to target precision `T`.

Parameters:
- `d`: Fourier beamlet propagator instance containing captured hits.
- `T`: Target floating-point number type.

Returns:
- Vector of Fourier k-points of type `Vector{FourierKPoint{T}}`.

Contracts & Invariants:
- Output vector elements are strictly type-stable with element type `FourierKPoint{T}`.
"""
function _extract_or_generate_kpoints(d::FourierBeamletPropagator, ::Type{T})::Vector{FourierKPoint{T}} where {T <: Real}
    if d isa FourierBeamletPropagator{T}
        return _extract_or_generate_kpoints(d)
    end
    pts_gen = generate_fourier_spectrum(d)
    return FourierKPoint{T}[FourierKPoint{T}(p.k4, p.E, p.weight) for p in pts_gen]
end

"""
    synthesize_field(d::FourierBeamletPropagator, grid::SpatialGrid{T, D}) where {T <: Real, D}

Synthesizes the complex vector electric field on the spatial detector grid from captured Fourier k-points.

Parameters:
- `d`: Propagator instance of type `FourierBeamletPropagator` containing detector hits or spectral k-points.
- `grid`: Spatial detector grid of type `SpatialGrid{T, D}` defining evaluation coordinates.

Returns:
- Complex electric field vector array of type `Array{SVector{3, Complex{T}}, D}` in V/m.

Contracts & Invariants:
- Array dimensions match the underlying spatial grid dimensions.
- Enforces electric field transversality and anti-padding invariant.
"""
function synthesize_field(d::FourierBeamletPropagator, grid::SpatialGrid{T, D})::Array{SVector{3, Complex{T}}, D} where {T <: Real, D}
    pts = _extract_or_generate_kpoints(d, T)
    dims = map(length, grid.ranges)
    E_out = zeros(SVector{3, Complex{T}}, dims)
    return _synthesize_field_nufft!(E_out, pts, grid, d)
end

"""
    synthesize_field(d::FourierBeamletPropagator{T}, ranges::NTuple{D, AbstractRange}) where {T <: Real, D}

Synthesizes the complex vector electric field on coordinates specified by range tuples.

Parameters:
- `d`: Propagator instance of type FourierBeamletPropagator{T} containing detector hits or spectral k-points.
- `ranges`: Coordinate ranges for each spatial dimension of type NTuple{D, AbstractRange}.

Output:
- Complex electric field vector array of type Array{SVector{3, Complex{T}}, D} in V/m.

Contracts & Invariants:
- All range lengths must be at least 4.
- Array dimensions match the lengths of the coordinate ranges.
"""
function synthesize_field(d::FourierBeamletPropagator{T}, ranges::NTuple{D, AbstractRange})::Array{SVector{3, Complex{T}}, D} where {T <: Real, D}
    return synthesize_field(d, SpatialGrid{T, D}(ranges))
end

export SpatialGrid, synthesize_field


