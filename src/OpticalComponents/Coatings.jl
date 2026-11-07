"""
    CoatingLayer(shape, model)

A free-standing thin layer, e.g. a dichroic film or a pellicle, that carries the
[`AbstractSurfaceModel`](@ref) `model` at every boundary at which it lies, see [`surface_model`](@ref).
The front of the layer is the side into which the normal of its `shape` points: the layers of a
[`ThinFilmCoating`](@ref) are listed from its back to its front. A coating on the face of a lens or a
prism is set with [`surface_model!`](@ref) instead.

```julia
pellicle = CoatingLayer(BeamletOptics.RectangularFlatMesh(25e-3, 25e-3), BeamsplitterCoating(0.08))
```
"""
struct CoatingLayer{T, S <: AbstractShape{T}, M <: AbstractSurfaceModel} <: AbstractObject{T}
    shape::S
    model::M
end

surface_model(layer::CoatingLayer) = layer.model

"""
    AbstractThinElement <: AbstractSurfaceModel

An infinitesimally thin element of a boundary, e.g. a lossless coating, a retarder film or a reflective
polarizer. It is described by its Jones matrices in reflection and transmission for light from both
of its sides, see [`ThinScattering`](@ref), without the Fresnel reflection of the bare boundary; the
transmitted field is scaled for the medium behind the boundary. Thin elements stack, see
[`SurfaceStack`](@ref).

# Implementation reqs.

## Fields

None required.

## Getters/setters

None required.

## Functions

- `scattering(element, frame::IncidenceFrame)`: the [`ThinScattering`](@ref) of the element for the
  ray of the [`IncidenceFrame`](@ref) `frame`
- `splitting_of(element)`: optional, see [`splitting_of`](@ref)

# Reflection phase

A lossless element reflects with the amplitudes ``r`` from its front and ``r' = -r`` from its back
(Stokes relations): only then do both outputs of an interferometer sum to its input. The front is the
side of [`isfacing`](@ref), see [`mirror_jones`](@ref) and
[Reflection phase at a beamsplitter coating](@ref).
"""
abstract type AbstractThinElement <: AbstractSurfaceModel end

"""
    IncidenceFrame{T, Tr}

The geometry of a ray at a thin element, from which an element computes its Jones matrices. They are
given in the s- and p-basis of the incident direction `d` and the reflected direction `r`; `s` is
perpendicular to the plane of incidence, `p_d = d × s` and `p_r = r × s`.

# Fields

- `d`, `r`: the incident and the reflected direction
- `s`, `p_d`, `p_r`: the s- and p-vectors of both directions
- `facing`: whether the ray comes from the front of the element, see [`isfacing`](@ref)
- `t`: the [`Transition`](@ref) of the boundary; its carrier gives the local frame of the axes of an
  element, see [`jones_angle`](@ref)
- `λ`: the wavelength in \\[m\\]
- `θ`: the angle of incidence in \\[rad\\], in the medium before the boundary
- `θ_air`: the angle of incidence in air of the same Snell invariant `n sin θ`, `π/2` beyond it
"""
struct IncidenceFrame{T, Tr <: Transition{T}}
    d::SVector{3, T}
    r::SVector{3, T}
    s::SVector{3, T}
    p_d::SVector{3, T}
    p_r::SVector{3, T}
    facing::Bool
    t::Tr
    λ::T
    θ::T
    θ_air::T
end

function IncidenceFrame(t::Transition, ray::AbstractRay{T}) where {T}
    d = SVector{3, T}(direction(ray)...)
    normal = normal3d(t)
    r = SVector{3, T}(reflection3d(d, normal)...)
    # The s-vector of the channels, see transfer
    s = s_vector(d, r, SVector{3, T}(normal3d(intersection(ray))...))
    θ = acos(min(abs(cos_angle3d(d, normal)), one(T)))
    θ_air = asin(min(refractive_index_in(t) * sin(θ), one(T)))
    return IncidenceFrame{T, typeof(t)}(d, r, s, cross(d, s), cross(r, s), isfacing(t, ray), t,
        T(wavelength(ray)), θ, θ_air)
end

scattering(element::AbstractThinElement, ::IncidenceFrame) =
    error(lazy"scattering not defined for $(typeof(element))")

"""
    ThinScattering{T}

The Jones matrices of a thin element for a ray, in the bases of the [`IncidenceFrame`](@ref): `r` and
`t` for light along the incident direction from the side of the ray, `r′` and `t′` for light along the
reflected direction from the other side.
"""
struct ThinScattering{T}
    r::SMatrix{2, 2, Complex{T}, 4}
    t::SMatrix{2, 2, Complex{T}, 4}
    r′::SMatrix{2, 2, Complex{T}, 4}
    t′::SMatrix{2, 2, Complex{T}, 4}
end

ThinScattering(r::SMatrix{2, 2, Complex{T}}, t::SMatrix{2, 2}, r′::SMatrix{2, 2},
    t′::SMatrix{2, 2}) where {T} = ThinScattering{T}(r, t, r′, t′)

"""
    cascade(near::ThinScattering, far::ThinScattering)

The scattering of two thin elements on top of each other, of which the ray meets `near` first. The
multiple reflections between them add coherently, without a phase between the elements (the Redheffer
star product).
"""
function cascade(near::ThinScattering, far::ThinScattering)
    # Round trips between the elements
    A = inv(I - near.r′ * far.r)
    B = inv(I - far.r * near.r′)
    return ThinScattering(near.r + near.t′ * far.r * A * near.t, far.t * A * near.t,
        far.r′ + far.t * near.r′ * B * far.t′, near.t′ * B * far.t′)
end

"""
    mirror_jones(frame::IncidenceFrame)

The Jones matrix of a lossless reflection with unit amplitude from the side of the ray of the `frame`:
`[1 0; 0 -1]` from the front, `[-1 0; 0 1]` from the back, like the external and the internal
reflection of the [`fresnel_coefficients`](@ref) at normal incidence.
"""
mirror_jones(frame::IncidenceFrame{T}) where {T} =
    SMatrix{2, 2, Complex{T}, 4}(frame.facing ? (1, 0, 0, -1) : (-1, 0, 0, 1))

"""
    jones_angle(frame, axis, p)

The angle of the `axis` of an element, in the local frame of its part, in the s- and p-basis of the
`frame` with the p-vector `p`: of its projection onto the plane transverse to the ray.
"""
function jones_angle(frame::IncidenceFrame, axis::AbstractVector, p::AbstractVector)
    a = orientation(part_of(carrier(frame.t))) * axis
    return atan(dot(a, p), dot(a, frame.s))
end

# The 2×2 matrix of the rotation by `φ`
rotation_jones(φ::T) where {T} = SMatrix{2, 2, T, 4}(cos(φ), sin(φ), -sin(φ), cos(φ))

# The projector onto the axis at the angle `φ`
projector_jones(φ) = (u = SVector(cos(φ), sin(φ)); u * transpose(u))

function channels(element::AbstractThinElement, t::Transition, ray::AbstractRay, cs)
    transmitted, reflected = cs
    # Total internal reflection: the ray does not reach the element
    islit(transmitted) || return cs
    frame = IncidenceFrame(t, ray)
    S = scattering(element, frame)
    # Transmitted field scaled for the medium behind the boundary
    normal = normal3d(t)
    cos_in = abs(cos_angle3d(direction(ray), normal))
    cos_out = abs(cos_angle3d(direction(transmitted), normal))
    f = sqrt(refractive_index_in(t) * cos_in / (refractive_index_out(t) * cos_out))
    a_t = channel_amplitude(ray, frame, S.t, jones_phase(S.t, 1))
    a_r = channel_amplitude(ray, frame, S.r, jones_phase(S.r, -1))
    return (transfer(transmitted, ray, jones_basis(f * S.t), a_t), transfer(reflected, ray, jones_basis(S.r), a_r))
end

# The 3×3 local Jones basis of the 2×2 matrix `J` in the s- and p-basis
jones_basis(J::SMatrix{2, 2}) =
    LocalJonesBasis(SMatrix{3, 3}(J[1, 1], J[2, 1], 0, J[1, 2], J[2, 2], 0, 0, 0, 1))

"""
    channel_amplitude(ray, frame, J, phase)

The amplitude factor of a channel with the Jones matrix `J` in the s- and p-basis of the `frame`, of
the unit phasor `phase`. Its magnitude is the square root of the power fraction of the channel: of
unpolarized light for a ray without a field, of the field of a [`PolarizedRay`](@ref) otherwise, such
that a channel that blocks the field carries no light. For a channel given by its power fractions
`Fs` and `Fp`, `J = diag(√Fs, √Fp)`.
"""
channel_amplitude(::AbstractRay, ::IncidenceFrame, J::AbstractMatrix, phase) =
    sqrt(sum(abs2, J) / 2) * phase

function channel_amplitude(ray::PolarizedRay, frame::IncidenceFrame, J::AbstractMatrix, phase)
    E = polarization(ray)
    E_sp = SVector(dot(frame.s, E), dot(frame.p_d, E))
    return norm(J * E_sp) / norm(E_sp) * phase
end

# The phase of the s- and p-amplitudes of `J`, in the sign convention `σ` of the s-component: -1 for a
# reflection, see mirror_jones
jones_phase(J::AbstractMatrix, σ) = phase_of(σ * J[1, 1] + J[2, 2])

# The matrix of the amplitudes of the power fractions `Fs` and `Fp`
power_jones(Fs, Fp) = SMatrix{2, 2}(sqrt(Fs), 0, 0, sqrt(Fp))

"""
    BeamsplitterCoating(R = 0.5)

A lossless, non-polarizing coating of the power reflectance `R` that splits a ray into a reflected
and a transmitted beam, e.g. `BeamsplitterCoating()` for 50:50 or `BeamsplitterCoating(0.7)` for
70:30. It acts between the media of the boundary, see [`Transition`](@ref): the transmitted beam is
refracted into the medium behind it, e.g. the substrate that it is coated on, the reflected beam
stays in the medium before it. Beyond the critical angle of these media, the ray is totally
reflected and not split. Its reflection phase follows the convention of
[`AbstractThinElement`](@ref).

A [`ThinBeamsplitter`](@ref) carries this coating, see [`surface_model`](@ref); see also
[`MirrorCoating`](@ref) and [`AntireflectionCoating`](@ref).

# Fields

- `reflectance`, `transmittance`: the field amplitudes `√R` and `√(1 - R)`
"""
struct BeamsplitterCoating{T <: Real} <: AbstractThinElement
    reflectance::T
    transmittance::T
end

function BeamsplitterCoating(R::Real = 0.5)
    0 ≤ R ≤ 1 ||
        throw(ArgumentError(lazy"the reflectance of a BeamsplitterCoating must lie in [0, 1], not $R"))
    r = sqrt(float(R))
    return BeamsplitterCoating(r, sqrt(1 - r^2))
end

reflectance(coating::BeamsplitterCoating) = coating.reflectance
transmittance(coating::BeamsplitterCoating) = coating.transmittance

splitting_of(::BeamsplitterCoating) = AllChannels()

scattering(coating::BeamsplitterCoating, frame::IncidenceFrame) =
    lossless_scattering(frame, reflectance(coating), reflectance(coating), transmittance(coating),
        transmittance(coating))

"""
    lossless_scattering(frame, rs, rp, τs, τp)

The [`ThinScattering`](@ref) of a lossless element that reflects the s- and p-components with the real
amplitudes `rs` and `rp` and transmits them with `τs` and `τp`, with the reflection phase of
[`mirror_jones`](@ref).
"""
function lossless_scattering(frame::IncidenceFrame{T}, rs, rp, τs, τp) where {T}
    M = mirror_jones(frame)
    ρ = SMatrix{2, 2, Complex{T}, 4}(rs, 0, 0, rp)
    τ = SMatrix{2, 2, Complex{T}, 4}(τs, 0, 0, τp)
    return ThinScattering(M * ρ, τ, -M * ρ, τ)
end

"""
    MirrorCoating()

An ideal, lossless mirror coating: it reflects all light, with the phase convention of a
[`BeamsplitterCoating`](@ref).
"""
MirrorCoating() = BeamsplitterCoating(1)

"""
    AntireflectionCoating()

An ideal anti-reflection coating: it transmits all light into the medium behind the boundary, without
the Fresnel reflection of the bare boundary.
"""
AntireflectionCoating() = BeamsplitterCoating(0)

"""
    SpectralCoating(R; splitting = AllChannels())
    SpectralCoating(Rs, Rp; splitting = AllChannels())
    SpectralCoating(λs, Rs; splitting = AllChannels())

A lossless coating with the power reflectance `R(λ, θ)` of the wavelength `λ` in \\[m\\] and the
angle of incidence in air `θ` in \\[rad\\], e.g. a dichroic mirror, an edge filter or a datasheet
curve; it transmits `1 - R`. A thin film responds to the Snell invariant `n sin θ` of the ray only,
the same on both of its sides, hence `θ` is the angle in air of this invariant, as datasheets give it:
for a ray at the angle `θ₁` in a medium of the index `n₁`, `sin θ = n₁ sin θ₁`, and `θ = π/2` for a
ray beyond the critical angle of air. With two functions, `Rs` and `Rp` are the reflectances of the s- and
p-components, e.g. of a polarizing beamsplitter coating. A table of the reflectances `Rs` at the
wavelengths `λs` is interpolated linearly and does not depend on the angle; beyond the table its
first or last value holds.

The phase convention is that of a [`BeamsplitterCoating`](@ref). With `splitting = Dominant()` only
the stronger part of the light is traced, see [`splitting_of`](@ref).

```julia
longpass = SpectralCoating((λ, θ) -> λ < 800e-9 ? 0.98 : 0.02)
```
"""
struct SpectralCoating{F, G, S <: AbstractSplitting} <: AbstractThinElement
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

function scattering(coating::SpectralCoating, frame::IncidenceFrame)
    Rs, Rp = coating.Rs(frame.λ, frame.θ_air), coating.Rp(frame.λ, frame.θ_air)
    return lossless_scattering(frame, sqrt(Rs), sqrt(Rp), sqrt(1 - Rs), sqrt(1 - Rp))
end

"""
    PolarizingCoating(Rs = 1, Rp = 0)

An ideal polarizing beamsplitter coating: a lossless [`SpectralCoating`](@ref) that reflects the
s-component of the field with the power reflectance `Rs` and the p-component with `Rp`, at all
wavelengths and angles. The s-component is perpendicular to the plane of incidence. Beams without
polarization are split with the mean reflectance.
"""
function PolarizingCoating(Rs::Real = 1, Rp::Real = 0)
    0 ≤ Rs ≤ 1 && 0 ≤ Rp ≤ 1 ||
        throw(ArgumentError(lazy"the reflectances of a PolarizingCoating must lie in [0, 1], not $Rs, $Rp"))
    return SpectralCoating((λ, θ) -> Rs, (λ, θ) -> Rp)
end

"""
    Retarder(Γ, fast_axis; λ0 = nothing)

An ideal retarder film, e.g. a quarter-wave (`Γ = π/2`) or half-wave (`Γ = π`) film: it transmits all
light and delays the component along the slow axis by the retardance `Γ` in \\[rad\\] behind the
component along the `fast_axis`, a vector in the local frame of the part that carries it, e.g.
`[1, 0, 0]` on a lens face whose normal is the local y-axis. At oblique incidence, the axis is projected
onto the plane transverse to the ray. With the design wavelength `λ0` in \\[m\\], the retardance is
`Γ λ0 / λ`, that of a zero-order plate without dispersion of its birefringence; otherwise it is the
same for all wavelengths.
"""
struct Retarder{T <: Real, L <: Nullable{Real}} <: AbstractThinElement
    Γ::T
    fast_axis::SVector{3, T}
    λ0::L
end

function Retarder(Γ::Real, fast_axis::AbstractVector{<:Real}; λ0::Nullable{Real} = nothing)
    T = float(promote_type(typeof(Γ), eltype(fast_axis)))
    return Retarder(T(Γ), SVector{3, T}(normalize(fast_axis)), isnothing(λ0) ? nothing : T(λ0))
end

"""
    fast_axis(layer::CoatingLayer)

Returns the unit vector, in global coordinates, of the fast axis of the [`Retarder`](@ref) of the
`layer`, e.g. of a [`QuarterWavePlate`](@ref). Its sign is arbitrary, since it is an axis.
"""
fast_axis(layer::CoatingLayer{T, S, <:Retarder}) where {T, S} =
    normalize(orientation(layer) * surface_model(layer).fast_axis)

"""Returns the retardance of the `retarder` in \\[rad\\] at the wavelength `λ` in \\[m\\]."""
retardance(retarder::Retarder{T, Nothing}, ::Real) where {T} = retarder.Γ
retardance(retarder::Retarder{T, <:Real}, λ::Real) where {T} = retarder.Γ * retarder.λ0 / λ

function scattering(retarder::Retarder, frame::IncidenceFrame{T}) where {T}
    Γ = retardance(retarder, frame.λ)
    J(φ) = rotation_jones(φ) * SMatrix{2, 2, Complex{T}, 4}(1, 0, 0, cis(Γ)) * rotation_jones(-φ)
    O = zero(SMatrix{2, 2, Complex{T}, 4})
    return ThinScattering(O, J(jones_angle(frame, retarder.fast_axis, frame.p_d)), O,
        J(jones_angle(frame, retarder.fast_axis, frame.p_r)))
end

"""
    ReflectivePolarizer(transmission_axis)

An ideal reflective polarizer film, e.g. a wire grid or a birefringent multilayer film: it transmits
the component of the field along the `transmission_axis`, a vector in the local frame of the part that
carries it, and reflects the component perpendicular to it, without loss. At oblique incidence, the
axis is projected onto the plane transverse to the ray. The reflection follows the phase convention of
[`AbstractThinElement`](@ref).
"""
struct ReflectivePolarizer{T <: Real} <: AbstractThinElement
    transmission_axis::SVector{3, T}
end

ReflectivePolarizer(axis::AbstractVector{<:Real}) =
    ReflectivePolarizer(SVector{3, float(eltype(axis))}(normalize(axis)))

splitting_of(::ReflectivePolarizer) = AllChannels()

function scattering(polarizer::ReflectivePolarizer, frame::IncidenceFrame{T}) where {T}
    M = mirror_jones(frame)
    P_d = projector_jones(jones_angle(frame, polarizer.transmission_axis, frame.p_d))
    P_r = projector_jones(jones_angle(frame, polarizer.transmission_axis, frame.p_r))
    C = Complex{T}
    return ThinScattering(M * (I - P_d), C.(P_d), -M * (I - P_r), C.(P_r))
end

"""
    SurfaceStack(elements::AbstractThinElement...; splitting = AllChannels())

Thin elements on top of each other at one boundary, e.g. a retarder film laminated under a reflective
polarizer, listed from the substrate outward: from the inside of a coated face of a body (see
[`surface_model!`](@ref)), or from the back of a [`CoatingLayer`](@ref) to its front. A ray meets them
in this order from the substrate and in the reverse order from outside; light reflected by an outer
element passes the inner ones again, and the multiple reflections between the elements add coherently,
see [`cascade`](@ref). A stack is itself a thin element and may contain stacks.

A [`ThinFilmCoating`](@ref) is no thin element: its layers form one coating, between the media of the
boundary.

```julia
qwp = Retarder(π / 2, [1, 0, 1])
fold = SurfaceStack(qwp, ReflectivePolarizer([1, 0, 0]))
```
"""
struct SurfaceStack{E <: Tuple, S <: AbstractSplitting} <: AbstractThinElement
    elements::E
    splitting::S
end

SurfaceStack(elements::AbstractThinElement...; splitting::AbstractSplitting = AllChannels()) =
    SurfaceStack(elements, splitting)

splitting_of(stack::SurfaceStack) = stack.splitting

function scattering(stack::SurfaceStack, frame::IncidenceFrame)
    # Elements in the order in which the ray meets them
    elements = frame.facing ? reverse(stack.elements) : stack.elements
    return foldl(cascade, map(e -> scattering(e, frame), elements))
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

A coherent stack of thin films, e.g. an anti-reflection, mirror, beamsplitter or filter coating of
known design. The `layers` are listed from the substrate outward: from the inside of a coated face
of a body (see [`surface_model!`](@ref)), or from the back of a [`CoatingLayer`](@ref) to its front.
The amplitudes of the reflected and transmitted s- and p-components follow from the characteristic
matrix of the stack (Born & Wolf, Principles of Optics, §1.6; Macleod, Thin-Film Optical Filters,
ch. 2), for both sides of incidence, any angle and absorbing layers; the media on both sides are
those of the boundary. Their phases follow from the stack as well. Without layers, the coating is
the bare boundary.

With `splitting = Dominant()` only the stronger part of the light is traced, e.g. for an
anti-reflection coating whose weak reflection is of no interest, see [`splitting_of`](@ref).

!!! info
    Beams without polarization, e.g. a [`GaussianBeamlet`](@ref), carry the power fractions of
    unpolarized light, with the phase of the sum of the s- and p-amplitudes; a polarized ray splits by
    the power fractions of its field.

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
    frame = IncidenceFrame(t, ray)
    a_r = channel_amplitude(ray, frame, power_jones(Rs, Rp), phase_of(rs + rp))
    a_t = channel_amplitude(ray, frame, power_jones(Ts, Tp), phase_of(ts + tp))
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
ghost images of a lens system. Beams without polarization carry the Fresnel power fractions of
unpolarized light. Each boundary splits the beam, hence a power cutoff (see
[`get_power_cutoff`](@ref)) or `depth_max` of [`solve_system!`](@ref) bounds the tree.
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
    frame = IncidenceFrame(t, ray)
    return (transfer(transmitted, transmitted.P, channel_amplitude(ray, frame, power_jones(Ts, Tp), 1)),
        transfer(reflected, reflected.P, channel_amplitude(ray, frame, power_jones(Rs, Rp), 1)))
end
