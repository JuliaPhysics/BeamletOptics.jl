"""
    CoatingLayer(shape, model)

A free-standing thin layer, e.g. a dichroic film or a pellicle, that carries the
[`AbstractSurfaceModel`](@ref) `model` at every boundary at which it lies, see [`surface_model`](@ref).
The front of the layer is the side into which the normal of its `shape` points: the layers of a
[`ThinFilmCoating`](@ref) are listed from its back to its front. A coating on the face of a lens or a
prism is set with [`surface_model!`](@ref) instead.

```julia
pellicle = CoatingLayer(BeamletOptics.RectangularFlatMesh(25e-3, 25e-3), BeamsplitterCoating(√0.08, √0.92))
```
"""
struct CoatingLayer{T, S <: AbstractShape{T}, M <: AbstractSurfaceModel} <: AbstractObject{T}
    shape::S
    model::M
end

surface_model(layer::CoatingLayer) = layer.model

"""
    BeamsplitterCoating{T} <: AbstractSurfaceModel

A lossless, non-polarizing coating that splits a ray into a transmitted and a reflected beam with the
field factors `reflectance` and `transmittance` (R² + T² = 1). It acts between the media of the
boundary, see [`Transition`](@ref): the transmitted beam is refracted into the medium behind it, e.g.
the substrate that it is coated on, the reflected beam stays in the medium before it. Beyond the
critical angle of these media, the ray is totally reflected and not split.

A [`ThinBeamsplitter`](@ref) carries this coating, see [`surface_model`](@ref); see also
[`MirrorCoating`](@ref) and [`AntireflectionCoating`](@ref).
"""
struct BeamsplitterCoating{T <: Real} <: AbstractSurfaceModel
    reflectance::T
    transmittance::T
end

reflectance(coating::BeamsplitterCoating) = coating.reflectance
transmittance(coating::BeamsplitterCoating) = coating.transmittance

splitting_of(::BeamsplitterCoating) = AllChannels()

"""
    channels(coating::BeamsplitterCoating, t::Transition, ray, cs)

The channels of the `coating`: the transmitted channel is refracted into the medium behind the
boundary with the field factor `transmittance(coating)`, the reflected one stays in the medium before
it with `reflectance(coating)`, see [`lossless_channels`](@ref).
"""
function channels(coating::BeamsplitterCoating, t::Transition, ray::AbstractRay, cs)
    r, τ = reflectance(coating), transmittance(coating)
    return lossless_channels(t, ray, cs, (r, r), (τ, τ))
end

"""
    lossless_channels(t::Transition, ray, cs, (rs, rp), (τs, τp))

The channels of a coating that reflects the s- and p-components of the field with the real amplitudes
`rs` and `rp` and transmits them with `τs` and `τp`, from the channels `cs` of the bare boundary
of the transition `t`. Beyond the critical angle the ray is totally reflected before it reaches the coating,
and `cs` is returned.

The field of a transmitted [`PolarizedRay`](@ref) is scaled for the medium behind the boundary: the
coating transmits a fraction of the power, and the power of a beam is proportional to `n cos(θ) |E|²`.
Beams without polarization carry the amplitudes of unpolarized light.

# Reflection phase jump

The reflection phase jump is θᵣ = π for a ray that comes from the front of the coating (see
[`isfacing`](@ref)), and 0 from the other side, as required by the Stokes relations (r' = -r) for a lossless
coating: only then do both outputs of an interferometer sum to its input. The s- and p-components are
reflected with `-r` and `r` from the other side, like an internal Fresnel reflection. The coating is
assumed to have an unambiguous normal, i.e. a 2D mesh. See also
[Reflection phase at a beamsplitter coating](@ref).
"""
function lossless_channels(t::Transition, ray::AbstractRay, cs, (rs, rp), (τs, τp))
    transmitted, reflected = cs
    # Total internal reflection: the ray does not reach the coating
    islit(transmitted) || return cs
    s = isfacing(t, ray) ? -1 : 1
    # Transmitted field scaled for the medium behind the boundary
    normal = normal3d(t)
    cos_in = abs(cos_angle3d(direction(ray), normal))
    cos_out = abs(cos_angle3d(direction(transmitted), normal))
    f = sqrt(refractive_index_in(t) * cos_in / (refractive_index_out(t) * cos_out))
    a_r = s * sqrt((rs^2 + rp^2) / 2)
    a_t = sqrt((τs^2 + τp^2) / 2)
    return (transfer(transmitted, ray, SPBasis(f * τs, 0, 0, f * τp), a_t),
        transfer(reflected, ray, SPBasis(-s * rs, 0, 0, s * rp), a_r))
end

"""
    MirrorCoating()

An ideal, lossless mirror coating: it reflects all light, with the phase convention of a
[`BeamsplitterCoating`](@ref).
"""
MirrorCoating() = BeamsplitterCoating(1.0, 0.0)

"""
    AntireflectionCoating()

An ideal anti-reflection coating: it transmits all light into the medium behind the boundary, without
the Fresnel reflection of the bare boundary.
"""
AntireflectionCoating() = BeamsplitterCoating(0.0, 1.0)

"""
    SpectralCoating(R; splitting = AllChannels())
    SpectralCoating(Rs, Rp; splitting = AllChannels())
    SpectralCoating(λs, Rs; splitting = AllChannels())

A lossless coating with the power reflectance `R(λ, θ)` of the wavelength `λ` in \\[m\\] and the angle of
incidence `θ` in \\[rad\\], e.g. a dichroic mirror, an edge filter or a datasheet curve; it transmits
`1 - R`. With two functions, `Rs` and `Rp` are the reflectances of the s- and p-components, e.g. of a
polarizing beamsplitter coating. A table of the reflectances `Rs` at the wavelengths `λs` is
interpolated linearly and does not depend on the angle; beyond the table its first or last value holds.

The phase convention is that of a [`BeamsplitterCoating`](@ref), see [`lossless_channels`](@ref). With
`splitting = Dominant()` only the stronger part of the light is traced, see [`splitting_of`](@ref).

```julia
longpass = SpectralCoating((λ, θ) -> λ < 800e-9 ? 0.98 : 0.02)
```
"""
struct SpectralCoating{F, G, S <: AbstractSplitting} <: AbstractSurfaceModel
    Rs::F
    Rp::G
    splitting::S
end

SpectralCoating(Rs, Rp; splitting::AbstractSplitting = AllChannels()) =
    SpectralCoating(Rs, Rp, splitting)
SpectralCoating(R; splitting::AbstractSplitting = AllChannels()) = SpectralCoating(R, R; splitting)
SpectralCoating(λs::AbstractVector{<:Real}, Rs::AbstractVector{<:Real}; kwargs...) =
    SpectralCoating(LinearTable(λs, Rs); kwargs...)

splitting_of(coating::SpectralCoating) = coating.splitting

function channels(coating::SpectralCoating, t::Transition, ray::AbstractRay, cs)
    λ = wavelength(ray)
    θ = acos(min(abs(cos_angle3d(direction(ray), normal3d(t))), 1))
    Rs, Rp = coating.Rs(λ, θ), coating.Rp(λ, θ)
    return lossless_channels(t, ray, cs, (sqrt(Rs), sqrt(Rp)), (sqrt(1 - Rs), sqrt(1 - Rp)))
end

"""
    LinearTable(xs, ys)

The values `ys` at the sorted points `xs`, interpolated linearly when called as `table(x)` or
`table(x, θ)`, which ignores `θ`. Beyond `xs`, the first or last value holds.
"""
struct LinearTable{X <: AbstractVector, Y <: AbstractVector}
    xs::X
    ys::Y
    function LinearTable(xs::X, ys::Y) where {X <: AbstractVector, Y <: AbstractVector}
        length(xs) == length(ys) || throw(ArgumentError("LinearTable needs as many values as points"))
        issorted(xs) || throw(ArgumentError("LinearTable needs sorted points"))
        return new{X, Y}(xs, ys)
    end
end

function (table::LinearTable)(x::Real)
    xs, ys = table.xs, table.ys
    x ≤ first(xs) && return first(ys)
    x ≥ last(xs) && return last(ys)
    i = searchsortedlast(xs, x)
    w = (x - xs[i]) / (xs[i + 1] - xs[i])
    return (1 - w) * ys[i] + w * ys[i + 1]
end

(table::LinearTable)(x::Real, ::Real) = table(x)

"""
    ThinFilmLayer(n, d)

A layer of a [`ThinFilmCoating`](@ref): the material of the [`AbstractRefractiveIndex`](@ref) `n`,
absorbing with its [`extinction_coefficient`](@ref), and the thickness `d` in \\[m\\]. See also
[`QuarterWaveLayer`](@ref).
"""
struct ThinFilmLayer{N <: AbstractRefractiveIndex, T <: Real}
    n::N
    d::T
    function ThinFilmLayer(n::N, d::T) where {N <: AbstractRefractiveIndex, T <: Real}
        d ≥ 0 || throw(ArgumentError("the thickness of a ThinFilmLayer must not be negative"))
        return new{N, T}(n, d)
    end
end

refractive_index(layer::ThinFilmLayer) = layer.n
thickness(layer::ThinFilmLayer) = layer.d

"""
    QuarterWaveLayer(n, λ0)

A [`ThinFilmLayer`](@ref) of the material `n` whose optical thickness is a quarter of the design
wavelength `λ0` in \\[m\\] at normal incidence: `d = λ0 / (4 n(λ0))`.
"""
QuarterWaveLayer(n::AbstractRefractiveIndex, λ0::Real) = ThinFilmLayer(n, λ0 / (4 * n(λ0)))

"""
    ThinFilmCoating(layers::ThinFilmLayer...; splitting = AllChannels())

A coherent stack of thin films, e.g. an anti-reflection, mirror, beamsplitter or filter coating of known
design. The `layers` are listed from the substrate outward: from the inside of a coated face of a body
(see [`surface_model!`](@ref)), or from the back of a [`CoatingLayer`](@ref) to its front. The amplitudes of the reflected and transmitted s- and
p-components follow from the characteristic matrix of the stack (Born & Wolf, Principles of Optics,
§1.6; Macleod, Thin-Film Optical Filters, ch. 2), for both sides of incidence, any angle and absorbing
layers; the media on both sides are those of the boundary. Their phases follow from the stack as well.
Without layers, the coating is the bare boundary.

With `splitting = Dominant()` only the stronger part of the light is traced, e.g. for an
anti-reflection coating whose weak reflection is of no interest, see [`splitting_of`](@ref).

!!! info
    Beams without polarization, e.g. a [`GaussianBeamlet`](@ref), carry the power fractions of
    unpolarized light, with the phase of the sum of the s- and p-amplitudes.

```julia
λ0 = 550e-9
mgf2 = ConstantRefractiveIndex(1.38)
ar = ThinFilmCoating(QuarterWaveLayer(mgf2, λ0))
```
"""
struct ThinFilmCoating{L <: Tuple, S <: AbstractSplitting} <: AbstractSurfaceModel
    layers::L
    splitting::S
end

ThinFilmCoating(layers::ThinFilmLayer...; splitting::AbstractSplitting = AllChannels()) =
    ThinFilmCoating(layers, splitting)

splitting_of(coating::ThinFilmCoating) = coating.splitting

function channels(coating::ThinFilmCoating, t::Transition, ray::AbstractRay, cs)
    transmitted, reflected = cs
    # Layers in the order in which the ray meets them
    layers = isfacing(t, ray) ? reverse(coating.layers) : coating.layers
    cos0 = abs(cos_angle3d(direction(ray), normal3d(t)))
    (rs, rp, ts, tp), (Rs, Rp, Ts, Tp) = thin_film_coefficients(layers, refractive_index_in(t),
        refractive_index_out(t), cos0, wavelength(ray))
    a_r = sqrt((Rs + Rp) / 2) * phase_of(rs + rp)
    a_t = sqrt((Ts + Tp) / 2) * phase_of(ts + tp)
    return (transfer(transmitted, ray, SPBasis(ts, 0, 0, tp), a_t),
        transfer(reflected, ray, SPBasis(-rs, 0, 0, rp), a_r))
end

# The unit phasor of `z`, one for zero
phase_of(z::Number) = iszero(z) ? one(z) : z / abs(z)

"""
    thin_film_coefficients(layers, n0, ns, cos0, λ)

Returns the field amplitudes `(rs, rp, ts, tp)` and the power fractions `(Rs, Rp, Ts, Tp)` of a stack of
the `layers`, met in this order, between the media of the real indices `n0` before and `ns` behind it,
for the cosine `cos0` of the angle of incidence and the wavelength `λ` in \\[m\\]. The amplitudes follow
the conventions of [`fresnel_coefficients`](@ref), with the field `∝ exp(i(kz - ωt))`: without layers
they are the Fresnel coefficients of the bare boundary.
"""
function thin_film_coefficients(layers::Tuple, n0::Real, ns::Real, cos0::Real, λ::Real)
    s2 = n0^2 * (1 - cos0^2)
    # N cos(θ) in a medium of the complex index N, decaying or propagating forward
    ncos(N) = (x = N^2 - s2; sqrt(complex(real(x), abs(imag(x)))))
    c0, cs = complex(n0 * cos0), ncos(complex(ns))
    # Tilted admittances of s and p, Macleod eq. 2.62
    ηs(_, c) = c
    ηp(N, c) = N^2 / c
    Ms, Mp = characteristic_matrices(layers, ncos, ηs, ηp, λ)
    rs, τs = stack_amplitudes(Ms, ηs(n0, c0), ηs(ns, cs))
    rp, τp = stack_amplitudes(Mp, ηp(n0, c0), ηp(ns, cs))
    # Tangential into total field of p, and power fractions, Macleod eq. 2.84
    tp = iszero(cs) ? zero(τp) : τp * cos0 * ns / cs
    Ts = real(ηs(ns, cs)) / real(ηs(n0, c0)) * abs2(τs)
    Tp = real(ηp(ns, cs)) / real(ηp(n0, c0)) * abs2(τp)
    return (rs, rp, τs, tp), (abs2(rs), abs2(rp), Ts, Tp)
end

# The characteristic matrices of the layers for s and p, with the field ∝ exp(i(kz - ωt)): the
# complex conjugate of Macleod eq. 2.97
function characteristic_matrices(layers::Tuple, ncos, ηs, ηp, λ::Real)
    C = typeof(complex(float(λ)))
    I2 = SMatrix{2, 2, C, 4}(1, 0, 0, 1)
    return foldl(layers; init = (I2, I2)) do (Ms, Mp), layer
        N = complex_index(refractive_index(layer), λ)
        c = ncos(N)
        δ = 2π * c * thickness(layer) / λ
        M(η) = SMatrix{2, 2, C, 4}(cos(δ), -im * η * sin(δ), -im * sin(δ) / η, cos(δ))
        (Ms * M(ηs(N, c)), Mp * M(ηp(N, c)))
    end
end

# The amplitudes of reflection and of the tangential field in transmission of the matrix `M` of a stack
# between the admittances `η0` and `ηs`, Macleod eq. 2.64
function stack_amplitudes(M::SMatrix{2, 2}, η0, ηs)
    B, Cm = M * SVector(one(ηs), ηs)
    return (η0 * B - Cm) / (η0 * B + Cm), 2η0 / (η0 * B + Cm)
end

"""
    FresnelGhosts()

The bare boundary, of which both the refracted and the reflected light are traced, e.g. to find the
ghost images of a lens system. Beams without polarization carry the Fresnel power fractions of unpolarized
light. Each boundary splits the beam, hence the keywords `power_cutoff` or `depth_max` of
[`solve_system!`](@ref) bound the tree.
"""
struct FresnelGhosts <: AbstractSurfaceModel end

splitting_of(::FresnelGhosts) = AllChannels()

function channels(::FresnelGhosts, t::Transition, ray::AbstractRay, cs)
    transmitted, reflected = cs
    # Total internal reflection: a single channel
    islit(transmitted) || return cs
    cos0 = abs(cos_angle3d(direction(ray), normal3d(t)))
    _, (Rs, Rp, Ts, Tp) = thin_film_coefficients((), refractive_index_in(t),
        refractive_index_out(t), cos0, wavelength(ray))
    return (transfer(transmitted, transmitted.P, sqrt((Ts + Tp) / 2)),
        transfer(reflected, reflected.P, sqrt((Rs + Rp) / 2)))
end
