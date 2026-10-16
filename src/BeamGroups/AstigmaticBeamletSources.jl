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
    dir_n = normalize(dir)
    # Use provided basis or fallback to automatic orthogonal basis
    e1 = isnothing(basis) ? normal3d(dir_n) : normalize(basis[1])
    e2 = isnothing(basis) ? normalize(cross(dir_n, e1)) : normalize(basis[2])

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    Δ = D / n_grid

    xs = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)
    ys = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)

    for x in xs
        for y in ys
            offset = x * e1 + y * e2
            # Support vector for the beamlet axes
            if randomize_axes
                base_s = normal3d(dir_n)
                ortho_s = cross(dir_n, base_s)
                phi = 2π * rand(rng)
                local_support = base_s * cos(phi) + ortho_s * sin(phi)
            else
                local_support = nothing
            end
            b = AstigmaticGaussianBeamlet(pos + offset, dir_n, λ, w0s; P0 = P0 / n_grid^2, support = local_support)
            push!(beams, b)
        end
    end
    # the beamlets overlap in the source plane, where their fields add
    _normalize_coherent_power!(beams, [position(first(rays(b.c))) for b in beams], w0s)
    # group orientation from the first grid axis, orthogonalized w.r.t. dir (sampling unchanged)
    e1_o = isnothing(basis) ? e1 : sampling_basis(dir_n, basis[1], T)
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, e1_o, T))
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
    if num_rays < num_rings * 20
        throw(ErrorException("No. of beamlets should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    if θ ≥ pi
        throw(ErrorException("Spherical source opening half-angle θ must be < π"))
    end

    # Calculate optimal sub-waist for smooth far-field overlap, from the angular spacing of the
    # rings, see θ_NA below
    Δθ = θ / (num_rings - 1)
    # We want beamlet divergence θ_div = λ / (π * w0s) to be approx overlap * Δθ
    w0s = λ / (π * overlap * Δθ)

    dir_n = normalize(dir)
    b1 = sampling_basis(dir_n, basis, T)
    b2 = normalize(cross(dir_n, b1))
    θ_NA = LinRange(0, θ, num_rings)

    # Initial ray count for power distribution, normalized for the coherent superposition below
    total_rays = num_rays
    P_sub = P0 / total_rays

    beams = Vector{AstigmaticGaussianBeamlet{T}}()

    # the field of a beamlet along `d` of the power P_sub with the polarization and phase of E0
    E0_of(d) = isnothing(E0) ? nothing : _jones_field(E0, d, electric_field(2P_sub / (π * w0s^2)))

    # Central beamlet
    push!(beams, AstigmaticGaussianBeamlet(pos, dir_n, λ, w0s; P0 = P_sub, E0 = E0_of(dir_n)))
    num_rays -= 1

    # Calculate circumference weights
    ndirs = [rotate3d(b2, step(θ_NA) * i) * dir_n for i in eachindex(θ_NA[2:end])]
    circm = norm.(ndirs .- dot.(ndirs, Ref(dir_n)) .* Ref(dir_n)) .* Ref(2π)
    # the outermost ring represents only the inner half of its band (trapezoidal rule), such
    # that every beamlet represents about the same solid angle
    circm[end] /= 2
    total = sum(circm)
    ds = total / num_rays

    n_rays = round.(Int, circm / ds)
    n_rays[end] += (num_rays - sum(n_rays))

    for (i, ndir) in enumerate(ndirs)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir_n, dphi)
        cdir = ndir
        for _ in 1:numEl
            # Polarization is set via the E0 kwarg. We could try to project it properly,
            # but default linear works for paraxial-ish cones.
            # Support vector for the beamlet axes
            if randomize_axes
                base_s = normal3d(cdir)
                ortho_s = cross(cdir, base_s)
                phi = 2π * rand(rng)
                local_support = base_s * cos(phi) + ortho_s * sin(phi)
            else
                local_support = nothing
            end
            push!(beams, AstigmaticGaussianBeamlet(pos, cdir, λ, w0s; P0 = P_sub, E0 = E0_of(cdir), support = local_support))
            cdir = RotMat * cdir
        end
    end

    # the beamlets overlap in the far field, where their fields add
    _normalize_coherent_power!(beams, [direction(first(rays(b.c))) for b in beams], λ / (π * w0s))
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
    if num_rays < num_rings * 20
        throw(ErrorException("No. of beamlets should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    if θ_x ≥ pi/2 || θ_y ≥ pi/2
        throw(ErrorException("Elliptical source opening half-angles θ_x and θ_y must be < π/2"))
    end

    # Calculate optimal sub-waist
    # We use the geometric mean of the angles to estimate the average ring spacing
    θ_avg = sqrt(θ_x * θ_y)
    # the spacing of the rings, see ρs below
    Δθ = θ_avg / (num_rings - 1)
    w0s = λ / (π * overlap * Δθ)

    dir_n = normalize(dir)
    b1 = sampling_basis(dir_n, basis, T)
    b2 = normalize(cross(dir_n, b1))
    
    # We use tangent space to define the elliptical rings
    tx = tan(θ_x)
    ty = tan(θ_y)

    # Initial ray count for power distribution, normalized for the coherent superposition below
    total_rays = num_rays
    P_sub = P0 / total_rays

    beams = Vector{AstigmaticGaussianBeamlet{T}}()

    # the field of a beamlet along `d` of the power P_sub with the polarization and phase of E0
    E0_of(d) = isnothing(E0) ? nothing : _jones_field(E0, d, electric_field(2P_sub / (π * w0s^2)))

    # Central beamlet
    push!(beams, AstigmaticGaussianBeamlet(pos, dir_n, λ, w0s; P0 = P_sub, E0 = E0_of(dir_n)))
    num_rays -= 1

    # Ring normalized radii, the center counted as the innermost of the num_rings rings
    ρs = LinRange(0, 1, num_rings)[2:end]

    # Calculate approximate circumference of each elliptical ring to distribute rays
    # Ramanujan's approximation for ellipse circumference: C ≈ π [3(a+b) - sqrt((3a+b)(a+3b))]
    circm = zeros(num_rings - 1)
    for (i, ρ) in enumerate(ρs)
        a = ρ * tx
        b = ρ * ty
        circm[i] = π * (3*(a+b) - sqrt((3a+b)*(a+3b)))
    end
    # the outermost ring represents only the inner half of its band (trapezoidal rule), such
    # that every beamlet represents about the same area of the tangent plane
    circm[end] /= 2
    
    total_circm = sum(circm)
    ds = total_circm / num_rays

    n_rays = round.(Int, circm / ds)
    n_rays[end] += (num_rays - sum(n_rays))

    for (i, ρ) in enumerate(ρs)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        
        a = ρ * tx
        b = ρ * ty
        
        # We use equal steps in the parametric angle t.
        # This ensures perfect symmetry and robust coverage.
        ts = LinRange(0, 2π, numEl + 1)[1:end-1]
        for t in ts
            # Parametric equation of ellipse in tangent plane
            u = a * cos(t)
            v = b * sin(t)
            
            # The 3D direction vector
            cdir = normalize(dir_n + u * b1 + v * b2)
            
            # Support vector for the beamlet axes
            if randomize_axes
                base_s = normal3d(cdir)
                ortho_s = cross(cdir, base_s)
                phi = 2π * rand(rng)
                local_support = base_s * cos(phi) + ortho_s * sin(phi)
            else
                local_support = nothing
            end
            push!(beams, AstigmaticGaussianBeamlet(pos, cdir, λ, w0s; P0 = P_sub, E0 = E0_of(cdir), support = local_support))
        end
    end

    # the beamlets overlap in the far field, where their fields add
    _normalize_coherent_power!(beams, [direction(first(rays(b.c))) for b in beams], λ / (π * w0s))
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, b1, T))
end

"""
    _normalize_coherent_power!(beams, vs, width)

Scales the fields of `beams`, beamlets of equal power that start in phase, such that their coherent
superposition carries the sum of their powers. Their profiles are Gaussians of the `width` (the
1/e² half width of the intensity) about the vectors `vs` in a common plane: the positions with the
waist in the source plane of collimated beamlets, or the directions with the far-field divergence
of beamlets from one point. Where the beamlets overlap, their fields add: two beamlets of the power
`P` whose vectors differ by `δ` contribute `P exp(-δ² / (2 width²))` to the power of the
superposition per pair, i.e. it carries `P Σᵢⱼ exp(-|vᵢ - vⱼ|² / (2 width²))` instead of `N P`.
This holds for all beamlets, also at the edge of the source, where fewer neighbors overlap.
"""
function _normalize_coherent_power!(beams, vs, width)
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
    _jones_field(E0, d, a)

The field of the amplitude `a` [V/m] along the direction `d` with the polarization and the phase of
the Jones vector `E0`: its component transverse to `d`, normalized. Throws an `ArgumentError` for
an `E0` along `d`.
"""
function _jones_field(E0, d, a)
    # `dot` conjugates its first argument, `d` is real
    p = E0 .- dot(d, E0) .* d
    n = norm(p)
    n > sqrt(eps(float(real(eltype(p))))) * norm(E0) || throw(ArgumentError("the polarization E0 is parallel to the direction $(d) of a beamlet"))
    return (a / n) .* p
end
