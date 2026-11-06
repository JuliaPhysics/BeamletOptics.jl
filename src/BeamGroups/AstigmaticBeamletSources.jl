"""
    CollimatedGaussianBeamletSource(pos, dir, D, λ, w0s; n_grid=20, P0)

Generates an `n_grid × n_grid` square array of `AstigmaticGaussianBeamlet`s, all pointing
in the same `dir`ection and having the same uniform amplitude. This is used to model
macroscopic plane waves or flat-top beams passing through hard apertures (e.g., for
Fraunhofer or Fresnel diffraction).

# Arguments
- `pos`: Center position of the source plane.
- `dir`: Propagation direction.
- `D`: Total width/height of the square aperture.
- `λ`: Wavelength.
- `w0s`: Sub-waist of each individual beamlet. For smooth overlap, `w0s ≈ D / n_grid`.
- `n_grid`: Number of beamlets along one axis (default `20` yields `400` total beamlets).
- `basis`: Optional tuple `(ex, ey)` of the macroscopic sampling grid axes, i.e. the 3D directions corresponding to the `x` and `y` grid axes. `ex` must not be zero or parallel to `dir`; its component normal to `dir` becomes the local x-axis of the group [`orientation`](@ref).
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical grid-alignment biases and is essential for preserving rotational symmetry in focused spots (e.g. Airy disks).
- `rng`: Random number generator to use for `randomize_axes`.
- `P0`: Total power of the source in [W]. The beamlets overlap and add coherently; their fields
  are scaled such that the superposition carries `P0`.
"""
function CollimatedGaussianBeamletSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        D::D2,
        λ::L,
        w0s::W;
        n_grid::Int = 20,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        P0::Real = get_default_power()
) where {P <: Real, D1 <: Real, D2 <: Real, L <: Real, W <: Real}
    T = promote_type(P, D1, D2, L, W)
    dir_n = normalize(T.(dir))
    e1, e2 = grid_axes(dir_n, basis)

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    Δ = D / n_grid

    xs = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)
    ys = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)

    for x in xs
        for y in ys
            offset = x * e1 + y * e2
            support = randomize_axes ? random_support(dir_n, rng) : nothing
            b = AstigmaticGaussianBeamlet(pos + offset, dir_n, λ, w0s; P0 = P0 / n_grid^2, support)
            push!(beams, b)
        end
    end
    # Normalize coherent power in the source plane
    normalize_coherent_power!(beams, [position(b) for b in beams], w0s)
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, sampling_basis(dir_n, basis, T), T))
end

"""
    SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings=10, num_rays=100*num_rings)

Decomposes a macroscopic spherical wave (or a highly divergent/focused beam) into a cone
of `AstigmaticGaussianBeamlet`s originating from a single point `pos`.

The angular spread of the beamlets is bounded by the half-angle `θ`. To ensure a smooth
far-field interference pattern without speckle, the sub-waist `w0s` of each beamlet is
automatically calculated such that their far-field divergence perfectly overlaps with
adjacent beamlets in the grid.

# Arguments
- `pos`: Origin point of the spherical wave.
- `dir`: Central propagation direction.
- `θ`: Half spread angle in radians.
- `λ`: Wavelength.
- `num_rings`: Number of concentric angular rings; the number of rings counts the center as the innermost one, i.e. `num_rings - 1` rings lie around it.
- `num_rays`: Total number of beamlets to generate.
- `overlap`: Scaling factor for the sub-waist divergence (default `1.2` ensures smooth overlap).
- `basis`: Optional reference vector (e.g. `[1,0,0]`) to define the starting azimuthal angle for the source rings.
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical artifacts in the far-field focus.
- `rng`: Random number generator to use for `randomize_axes`.
- `P0`: Total power of the spherical source in [W]. The beamlets overlap and add coherently in the
  far field; their fields are scaled such that the superposition carries `P0`.
- `E0`: Optional Jones vector defining the polarization and phase of the spherical wave; its
  component transverse to each beamlet is used, its magnitude does not matter, the power is `P0`.

!!! info "Reproducible sampling"
    If no `basis` is passed, the orthogonal basis vectors are derived from `dir` deterministically,
    so two sources sharing the same arguments sample exactly the same beamlet directions. Pass a
    `basis` to rotate the azimuthal sampling of a source about its own axis. A `basis` that is
    zero or parallel to `dir` throws.
"""
function SphericalGaussianBeamletSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ::H,
        λ::L;
        num_rings::Int = 10,
        num_rays::Int = 100 * num_rings,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, AbstractVector} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        P0::Real = get_default_power(),
        E0 = nothing
) where {P <: Real, D <: Real, H <: Real, L <: Real}
    T = promote_type(P, D, H, L)
    num_rays ≥ 20num_rings ||
        throw(ArgumentError("SphericalGaussianBeamletSource needs num_rays ≥ 20 num_rings = $(20num_rings), got $num_rays"))
    θ < π || throw(ArgumentError("SphericalGaussianBeamletSource needs a half spread angle θ < π, got $θ"))

    # Sub-waist for smooth far-field overlap of adjacent rings
    Δθ = θ / (num_rings - 1)
    w0s = λ / (π * overlap * Δθ)

    dir_n = normalize(T.(dir))
    b1 = sampling_basis(dir_n, basis, T)
    dirs = cone_ring_directions(dir_n, b1, θ, num_rings, num_rays - 1)
    beams = cone_beamlets(pos, dir_n, dirs, λ, w0s, P0, E0, randomize_axes, rng, T)
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, b1, T))
end

"""
    EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ; num_rings=10, num_rays=100*num_rings, ...)

Spawns a coherent array of `AstigmaticGaussianBeamlet`s distributed on an elliptical cap.
This is ideal for modelling sources with different divergences in the fast and slow axes,
such as tapered amplifiers or edge-emitting laser diodes.

# Arguments
- `pos`: Origin point of the spherical wave.
- `dir`: Central propagation direction.
- `θ_x`, `θ_y`: Half spread angles in radians for the two principal axes.
- `λ`: Wavelength.
- `num_rings`: Number of concentric angular rings; the number of rings counts the center as the innermost one, i.e. `num_rings - 1` rings lie around it.
- `num_rays`: Total number of beamlets to generate.
- `overlap`: Scaling factor for the sub-waist divergence (default `1.2` ensures smooth overlap).
- `basis`: Optional reference vector (e.g. `[1,0,0]`) to define the starting azimuthal angle for the source rings.
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated.
- `rng`: Random number generator to use for `randomize_axes`.
- `P0`: Total power of the source in [W]. The beamlets overlap and add coherently in the far field;
  their fields are scaled such that the superposition carries `P0`.
- `E0`: Optional Jones vector defining the polarization and phase of the wave; its component
  transverse to each beamlet is used, its magnitude does not matter, the power is `P0`.

!!! info "Reproducible sampling"
    If no `basis` is passed, the orthogonal basis vectors are derived from `dir` deterministically,
    so two sources sharing the same arguments sample exactly the same beamlet directions. Pass a
    `basis` to rotate the azimuthal sampling of a source about its own axis. A `basis` that is
    zero or parallel to `dir` throws.
"""
function EllipticalGaussianBeamletSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ_x::H1,
        θ_y::H2,
        λ::L;
        num_rings::Int = 10,
        num_rays::Int = 100 * num_rings,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, AbstractVector} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        P0::Real = get_default_power(),
        E0 = nothing
) where {P <: Real, D <: Real, H1 <: Real, H2 <: Real, L <: Real}
    T = promote_type(P, D, H1, H2, L)
    num_rays ≥ 20num_rings ||
        throw(ArgumentError("EllipticalGaussianBeamletSource needs num_rays ≥ 20 num_rings = $(20num_rings), got $num_rays"))
    θ_x < π / 2 && θ_y < π / 2 ||
        throw(ArgumentError("EllipticalGaussianBeamletSource needs half spread angles θ_x, θ_y < π/2, got $θ_x, $θ_y"))

    # Sub-waist from the geometric mean of the ring spacings
    θ_avg = sqrt(θ_x * θ_y)
    Δθ = θ_avg / (num_rings - 1)
    w0s = λ / (π * overlap * Δθ)

    dir_n = normalize(T.(dir))
    b1 = sampling_basis(dir_n, basis, T)
    b2 = normalize(cross(dir_n, b1))

    # Elliptical rings in the tangent plane
    tx = tan(θ_x)
    ty = tan(θ_y)
    ρs = LinRange(0, 1, num_rings)[2:end]

    # Ring circumferences (Ramanujan)
    circm = zeros(num_rings - 1)
    for (i, ρ) in enumerate(ρs)
        a = ρ * tx
        b = ρ * ty
        circm[i] = π * (3*(a+b) - sqrt((3a+b)*(a+3b)))
    end
    counts = rays_per_ring(circm, num_rays - 1)

    # Equal steps in the parametric angle
    dirs = [normalize(dir_n + ρ * tx * cos(t) * b1 + ρ * ty * sin(t) * b2)
            for (ρ, n) in zip(ρs, counts) if !iszero(n) for t in LinRange(0, 2π, n + 1)[1:end-1]]
    beams = cone_beamlets(pos, dir_n, dirs, λ, w0s, P0, E0, randomize_axes, rng, T)
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, b1, T))
end

"""
    cone_beamlets(pos, dir, dirs, λ, w0s, P0, E0, randomize_axes, rng, T)

The beamlets of a cone source at `pos` of the total power `P0`: a central beamlet along `dir`
and one along each of `dirs`, of the waist `w0s`, normalized for their coherent superposition in
the far field.
"""
function cone_beamlets(
        pos::AbstractArray,
        dir::AbstractArray,
        dirs::AbstractVector,
        λ::Real,
        w0s::Real,
        P0::Real,
        E0,
        randomize_axes::Bool,
        rng::AbstractRNG,
        ::Type{T}
) where {T <: Real}
    P_sub = P0 / (length(dirs) + 1)
    a = electric_field(2P_sub / (π * w0s^2))
    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    # Central beamlet
    push!(beams, AstigmaticGaussianBeamlet(pos, dir, λ, w0s; P0 = P_sub, E0 = jones_field(E0, dir, a)))
    for d in dirs
        support = randomize_axes ? random_support(d, rng) : nothing
        push!(beams, AstigmaticGaussianBeamlet(pos, d, λ, w0s; P0 = P_sub, E0 = jones_field(E0, d, a), support))
    end
    # Normalize coherent power in the far field
    normalize_coherent_power!(beams, [direction(b) for b in beams], λ / (π * w0s))
    return beams
end

"""
    normalize_coherent_power!(beams, vs, width)

Scales the fields of `beams`, of equal power and in phase, such that their coherent superposition
carries the sum of their powers; `vs` are their positions or directions, `width` the 1/e² half
width of their profiles.
"""
function normalize_coherent_power!(beams::AbstractVector{<:AstigmaticGaussianBeamlet}, vs, width)
    # Truncate overlap at 8 widths
    cutoff = (8width)^2
    S = zero(width)
    for i in eachindex(vs), j in eachindex(vs)
        δ2 = sum(abs2, vs[i] - vs[j])
        δ2 < cutoff && (S += exp(-δ2 / (2width^2)))
    end
    scale = sqrt(length(vs) / S)
    for b in beams
        electric_field!(b, electric_field(b) * scale)
    end
    return beams
end

"""
    jones_field(E0, d, a)

The field of the amplitude `a` [V/m] along the direction `d` with the polarization and the phase of
the Jones vector `E0`: its component transverse to `d`, normalized. Throws an `ArgumentError` for
an `E0` along `d`; returns `nothing` for `E0 = nothing`.
"""
function jones_field(E0::AbstractVector, d, a)
    p = _transverse(E0, d)
    n = norm(p)
    # Parallel to d within round-off
    n > sqrt(eps(float(real(eltype(p))))) * norm(E0) || throw(ArgumentError("the polarization E0 is parallel to the direction $(d) of a beamlet"))
    return (a / n) .* p
end

jones_field(::Nothing, d, a) = nothing

"""
    grid_axes(dir, basis) -> (e1, e2)

The unit axes of a beamlet grid normal to the unit vector `dir`: `normal3d(dir)` and
`dir × e1` for `basis = nothing`, else the normalized axes of the tuple `basis`.
"""
function grid_axes(dir::AbstractVector, ::Nothing)
    e1 = normal3d(dir)
    return e1, normal3d(dir, e1)
end

grid_axes(::AbstractVector, basis::Tuple{AbstractVector, AbstractVector}) =
    (normalize(basis[1]), normalize(basis[2]))

"""
    random_support(dir, rng)

A unit vector normal to the unit vector `dir` at a random azimuth, from one draw of `rng`: the
support of the principal axes of a beamlet with randomized axes.
"""
function random_support(dir::AbstractVector, rng::AbstractRNG)
    base = normal3d(dir)
    ortho = cross(dir, base)
    φ = 2π * rand(rng)
    return base * cos(φ) + ortho * sin(φ)
end
