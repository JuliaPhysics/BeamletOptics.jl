"""
    AstigmaticGaussianBeamlet{T} <: AbstractBeam{T, PolarizedRay{T}}

Complex ray representation of a general **astigmatic** Gaussian beam as per the formalism described by A. Greynolds (1986) and N. Worku (2017).
The beam is described by a chief [`PolarizedRay`](@ref) and 8 auxiliary [`Ray`](@ref)s that encode the waist radius and divergence in two orthogonal transverse axes.
The beam quality `M2` is considered via the divergence angle. All equations are based on the following publications:

**Alan Greynolds, "Vector formulation of the ray-equivalent method for general Gaussian beam propagation." Curr. Dev. Opt. Eng. Diffraction Phenom. Vol. 679. SPIE, 1986**

and

**Norman Worku and Herbert Gross, "Vectorial field propagation through high NA objectives using polarized Gaussian beam decomposition." OTOM XIV. Vol. 10347. SPIE, 2017**

# Fields

- `c`:   chief [`Beam`](@ref) of [`PolarizedRay`](@ref)s
- `wxp`: waist  x-positive [`Beam`](@ref) of [`Ray`](@ref)s
- `wxm`: waist  x-negative [`Beam`](@ref) of [`Ray`](@ref)s
- `wyp`: waist  y-positive [`Beam`](@ref) of [`Ray`](@ref)s
- `wym`: waist  y-negative [`Beam`](@ref) of [`Ray`](@ref)s
- `dxp`: divergence x-positive [`Beam`](@ref) of [`Ray`](@ref)s
- `dxm`: divergence x-negative [`Beam`](@ref) of [`Ray`](@ref)s
- `dyp`: divergence y-positive [`Beam`](@ref) of [`Ray`](@ref)s
- `dym`: divergence y-negative [`Beam`](@ref) of [`Ray`](@ref)s
- `parent`:   reference to the parent beam (or `nothing`)
- `children`: vector of child beams (e.g. after beam-splitting)

# Additional information

!!! info "Waist and field calculation"
    For a given beamlet the Gaussian beam parameters can be obtained via the [`gauss_parameters`](@ref) function.
"""
mutable struct AstigmaticGaussianBeamlet{T, RC <: PolarizedRay{T}, RA <: Ray{T}} <: AbstractBeam{T, RC}
    c::Beam{T, RC}       # chief
    wxp::Beam{T, RA}     # waist x +
    wxm::Beam{T, RA}     # waist x -
    wyp::Beam{T, RA}     # waist y +
    wym::Beam{T, RA}     # waist y -
    dxp::Beam{T, RA}     # div. x +
    dxm::Beam{T, RA}     # div. x -
    dyp::Beam{T, RA}     # div. y +
    dym::Beam{T, RA}     # div. y -
    parent::Nullable{AstigmaticGaussianBeamlet{T, RC, RA}}
    children::Vector{AstigmaticGaussianBeamlet{T, RC, RA}}
 end

function AstigmaticGaussianBeamlet(
        c::Beam{T, RC},
        wxp::Beam{T, RA},
        wxm::Beam{T, RA},
        wyp::Beam{T, RA},
        wym::Beam{T, RA},
        dxp::Beam{T, RA},
        dxm::Beam{T, RA},
        dyp::Beam{T, RA},
        dym::Beam{T, RA}
) where {T <: Real, RC <: PolarizedRay{T}, RA <: Ray{T}}
    return AstigmaticGaussianBeamlet{T, RC, RA}(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym,
        nothing,
        Vector{AstigmaticGaussianBeamlet{T, RC, RA}}())
end

"""
    AstigmaticGaussianBeamlet(position, direction, λ, w0; kwargs...)

Constructs an astigmatic Gaussian beamlet at its waist with the specified beam parameters.
In the 4-argument version, the beam has circular symmetry with waist radius `w0`.
In the 5-argument version, independent waists `w0_x` and `w0_y` can be specified.

# Arguments

## Inputs

- `position`: origin of the beamlet
- `direction`: direction of the beamlet
- `λ`: wavelength of the beamlet in [m]. Default value is 1000 nm.
- `w0`: beam waist (radius) in [m]. Default value is 1 mm.
- `w0_x`, `w0_y`: independent beam waists in [m].

## Keyword Arguments

- `M2`, `M2_x`, `M2_y`: beam quality factors. Default is 1
- `P0`: beam total power in [W]. Default is 1 mW.
- `E0`: electric field vector in [V/m]. Default is `nothing` (aligned with support axes, scaled by `P0`).
- `support`: optional support vector for basis construction
- `z0`: beam waist offset in [m]. Default is 0 m

# Additional information

!!! warning "Optical invariant check"
    When using the `solve_system!` function, the beamlet invariant will be checked for each interaction. If the invariant is violated, tracing will be stopped.
"""
function AstigmaticGaussianBeamlet(
        position::AbstractArray,
        direction::AbstractArray,
        λ::Real = get_default_wavelength(),
        w0::Real = get_default_waist();
        # kwargs
        M2::Real = 1,
        P0::Real = get_default_power(),
        E0::Union{Nothing, AbstractArray{<:Number}} = nothing,
        support::Union{Nothing, AbstractArray{<:Real}} = nothing,
        z0::Real = 0
)
    return AstigmaticGaussianBeamlet(position, direction, λ, w0, w0;
        M2_x = M2, M2_y = M2, P0, E0, support, z0)
end

"""
    AstigmaticGaussianBeamlet(position, direction, λ, w0_x, w0_y; kwargs...)

Constructs an astigmatic Gaussian beamlet with independent waists `w0_x` and `w0_y`.
This constructor supports modeling astigmatic sources where the waists in X and Y are at different locations along the beam axis.

# Arguments
- `position`: Global reference point of the beamlet.
- `direction`: Propagation direction vector.
- `λ`: Wavelength.
- `w0_x`: Waist radius in the X direction (defined by the support vector).
- `w0_y`: Waist radius in the Y direction.

# Keyword Arguments
- `z0_x`: Position of the X waist relative to `position` along the propagation axis. Defaults to `z0`.
- `z0_y`: Position of the Y waist relative to `position` along the propagation axis. Defaults to `z0`.
- `z0`: Default waist position if `z0_x` and `z0_y` are not specified. Also defines the starting point of the chief ray as `position + z0 * direction`.
- `M2_x`, `M2_y`: Beam quality factors. Default is 1.
- `P0`: Total power in [W].
- `E0`: Optional Jones vector for polarization.
- `support`: Optional vector orthogonal to `direction` to define the X axis.

# Additional information

!!! warning "Optical invariant check"
    When using the `solve_system!` function, the beamlet invariant will be checked for each interaction. If the invariant is violated, tracing will be stopped.
"""
function AstigmaticGaussianBeamlet(
        position::AbstractArray,
        direction::AbstractArray,
        λ::Real,
        w0_x::Real,
        w0_y::Real;
        # kwargs
        M2_x::Real = 1,
        M2_y::Real = 1,
        P0::Real = get_default_power(),
        E0::Union{Nothing, AbstractArray{<:Number}} = nothing,
        support::Union{Nothing, AbstractArray{<:Real}} = nothing,
        z0::Real = 0,
        z0_x::Real = z0,
        z0_y::Real = z0
)
    # Create orthogonal vectors for construction purposes (right-handed)
    direction = normalize(direction)
    if isnothing(support)
        s1 = normal3d(direction)
    else
        if !isorthogonal3d(direction, support)
            error("Ray direction and support vector must be orthogonal!")
        end
        s1 = normalize(support)
    end
    s2 = cross(direction, s1)

    # Handle polarization and power
    if isnothing(E0)
        # Calculate E0 magnitude based on P0
        I0 = (2 * P0) / (π * w0_x * w0_y)
        E_phasor = electric_field(I0)
        E0 = s1 .* E_phasor # Default to linear polarization along first transverse axis
    else
        if !isorthogonal3d(direction, E0)
            # If user provided E0 but it's not orthogonal, project it.
            # This is a convenience for tilted setups.
            E0 = E0 .- dot(E0, direction) .* direction
        end
    end

    # Divergence angles
    θx = divergence_angle(λ, w0_x, M2_x)
    θy = divergence_angle(λ, w0_y, M2_y)
    # Waist rays
    wxp = Ray(position + s1 * w0_x + z0_x * direction, direction, λ)
    wxm = Ray(position - s1 * w0_x + z0_x * direction, direction, λ)
    wyp = Ray(position + s2 * w0_y + z0_y * direction, direction, λ)
    wym = Ray(position - s2 * w0_y + z0_y * direction, direction, λ)
    # Divergence rays
    div_dir_xp = normalize(direction + s1 * tan(θx))
    div_dir_xm = normalize(direction - s1 * tan(θx))
    div_dir_yp = normalize(direction + s2 * tan(θy))
    div_dir_ym = normalize(direction - s2 * tan(θy))
    # Corrected divergence ray positions to ensure waist is at z0_x and z0_y
    dxp = Ray(position - z0_x * s1 * tan(θx), div_dir_xp, λ)
    dxm = Ray(position + z0_x * s1 * tan(θx), div_dir_xm, λ)
    dyp = Ray(position - z0_y * s2 * tan(θy), div_dir_yp, λ)
    dym = Ray(position + z0_y * s2 * tan(θy), div_dir_ym, λ)
    # Chief ray
    c = PolarizedRay(position + z0 * direction, direction, λ, E0)
    return AstigmaticGaussianBeamlet(
        Beam(c),
        Beam(wxp),
        Beam(wxm),
        Beam(wyp),
        Beam(wym),
        Beam(dxp),
        Beam(dxm),
        Beam(dyp),
        Beam(dym)
    )
end

"""Return a tuple of all 9 component beams of the [`AstigmaticGaussianBeamlet`](@ref)."""
@inline _component_beams(agb::AstigmaticGaussianBeamlet) = (
    agb.c, agb.wxp, agb.wxm, agb.wyp, agb.wym, agb.dxp, agb.dxm, agb.dyp, agb.dym)

"""Return a tuple of the 8 auxiliary (parabasal) beams."""
@inline _aux_beams(agb::AstigmaticGaussianBeamlet) = (
    agb.wxp, agb.wxm, agb.wyp, agb.wym, agb.dxp, agb.dxm, agb.dyp, agb.dym)

# Implementation of AbstractBeam interface
rays(agb::AstigmaticGaussianBeamlet) = rays(agb.c)

Base.length(agb::AstigmaticGaussianBeamlet) = length(agb.c)
optical_path_length(agb::AstigmaticGaussianBeamlet) = optical_path_length(agb.c)

isentering(agb::AstigmaticGaussianBeamlet, id::Int) = isentering(rays(agb.c)[id])

function _reset_beam!(agb::AstigmaticGaussianBeamlet)
    for beam in _component_beams(agb)
        _reset_beam!(beam)
    end
    _drop_beams!(agb)
    return nothing
end

point_on_beam(agb::AstigmaticGaussianBeamlet, t::Real) = point_on_beam(agb.c, t)

wavelength(agb::AstigmaticGaussianBeamlet) = wavelength(first(rays(agb.c)))
direction(agb::AstigmaticGaussianBeamlet) = direction(first(rays(agb.c)))
position(agb::AstigmaticGaussianBeamlet) = position(first(rays(agb.c)))
polarization(agb::AstigmaticGaussianBeamlet) = polarization(first(rays(agb.c)))

electric_field(agb::AstigmaticGaussianBeamlet) = polarization(agb)
function electric_field!(agb::AstigmaticGaussianBeamlet, new_E0)
    polarization!(first(rays(agb.c)), new_E0)
end

function refractive_index(agb::AstigmaticGaussianBeamlet, id::Int)
    return refractive_index(rays(agb.c)[id])
end

function refractive_index!(agb::AstigmaticGaussianBeamlet, id::Int, n_new::Real)
    for beam in _component_beams(agb)
        refractive_index!(rays(beam)[id], n_new)
    end
    return nothing
end

"""
    parent!(child::AstigmaticGaussianBeamlet, parent::AstigmaticGaussianBeamlet)

Links parent for tree navigation and ensures the chief beam parent is also linked.
"""
function parent!(child::AstigmaticGaussianBeamlet, parent::AstigmaticGaussianBeamlet)
    child.parent = parent
    child.c.parent = parent.c
    return nothing
end

function Base.show(io::IO, agb::AstigmaticGaussianBeamlet)
    p0 = position(first(rays(agb.c)))
    d0 = direction(first(rays(agb.c)))
    λ = wavelength(agb)
    _, w1, w2 = waist_parameters(agb, 0.0)
    print(io,
        "AstigmaticGaussianBeamlet(pos: $p0, dir: $d0, λ: $λ, w_x: $(norm(w1)), w_y: $(norm(w2)))")
end

# AbstractTrees integration
AbstractTrees.children(agb::AstigmaticGaussianBeamlet) = agb.children
AbstractTrees.nodevalue(agb::AstigmaticGaussianBeamlet) = agb

"""
    AstigmaticGaussianBeamletInteraction <: AbstractInteraction

Stores the interaction result for all 9 component beams of an
[`AstigmaticGaussianBeamlet`](@ref). Uses the hint from the chief beam interaction.

# Fields

- `chief`: [`BeamInteraction`](@ref) for the chief ray
- `wxp`, `wxm`, `wyp`, `wym`: waist beam interactions
- `dxp`, `dxm`, `dyp`, `dym`: divergence beam interactions
"""
struct AstigmaticGaussianBeamletInteraction{R <: Real} <: AbstractInteraction
    chief::BeamInteraction{R}
    wxp::BeamInteraction{R}
    wxm::BeamInteraction{R}
    wyp::BeamInteraction{R}
    wym::BeamInteraction{R}
    dxp::BeamInteraction{R}
    dxm::BeamInteraction{R}
    dyp::BeamInteraction{R}
    dym::BeamInteraction{R}
end

hint(i::AstigmaticGaussianBeamletInteraction) = hint(i.chief)
function hint!(i::AstigmaticGaussianBeamletInteraction, new_hint::Nullable{Hint})
    hint!(i.chief, new_hint)
end

"""
    interact3d(system, object, agb::AstigmaticGaussianBeamlet, ray_id)

Generic dispatch: traces each component beam's ray through `object` independently.
Returns `nothing` if any interaction fails; otherwise returns an
[`AstigmaticGaussianBeamletInteraction`](@ref).
"""
function interact3d(system::AbstractSystem,
        object::AbstractObject,
        agb::AstigmaticGaussianBeamlet{R},
        ray_id::Int) where {R}
    i_c = interact3d(system, object, agb.c, rays(agb.c)[ray_id])
    i_wxp = interact3d(system, object, agb.wxp, rays(agb.wxp)[ray_id])
    i_wxm = interact3d(system, object, agb.wxm, rays(agb.wxm)[ray_id])
    i_wyp = interact3d(system, object, agb.wyp, rays(agb.wyp)[ray_id])
    i_wym = interact3d(system, object, agb.wym, rays(agb.wym)[ray_id])
    i_dxp = interact3d(system, object, agb.dxp, rays(agb.dxp)[ray_id])
    i_dxm = interact3d(system, object, agb.dxm, rays(agb.dxm)[ray_id])
    i_dyp = interact3d(system, object, agb.dyp, rays(agb.dyp)[ray_id])
    i_dym = interact3d(system, object, agb.dym, rays(agb.dym)[ray_id])
    if any(isnothing, (i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym))
        return nothing
    end
    return AstigmaticGaussianBeamletInteraction{R}(
        i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym)
end

function interact3d(system::AbstractSystem,
        mi::MultiIntersection,
        agb::AstigmaticGaussianBeamlet{R},
        ray_id::Int) where {R}
    c_obj = coating(mi)
    obj_en = entering(mi)
    obj_ex = exiting(mi)
    bs_obj = c_obj isa AbstractBeamsplitter ? c_obj : (obj_en isa AbstractBeamsplitter ? obj_en : (obj_ex isa AbstractBeamsplitter ? obj_ex : nothing))
    if bs_obj !== nothing
        ambient = ambient_medium(system)
        trans = resolve_transition(mi, ambient, rays(agb.c)[ray_id])
        return interact3d(system, bs_obj, trans, mi.hit, agb, ray_id)
    end

    i_c = interact3d(system, mi, agb.c, rays(agb.c)[ray_id])
    i_wxp = interact3d(system, mi, agb.wxp, rays(agb.wxp)[ray_id])
    i_wxm = interact3d(system, mi, agb.wxm, rays(agb.wxm)[ray_id])
    i_wyp = interact3d(system, mi, agb.wyp, rays(agb.wyp)[ray_id])
    i_wym = interact3d(system, mi, agb.wym, rays(agb.wym)[ray_id])
    i_dxp = interact3d(system, mi, agb.dxp, rays(agb.dxp)[ray_id])
    i_dxm = interact3d(system, mi, agb.dxm, rays(agb.dxm)[ray_id])
    i_dyp = interact3d(system, mi, agb.dyp, rays(agb.dyp)[ray_id])
    i_dym = interact3d(system, mi, agb.dym, rays(agb.dym)[ray_id])
    if any(isnothing, (i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym))
        return nothing
    end
    return AstigmaticGaussianBeamletInteraction{R}(
        i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym)
end

function Base.push!(agb::AstigmaticGaussianBeamlet{T},
        interaction::AstigmaticGaussianBeamletInteraction{T}) where {T}
    push!(agb.c, interaction.chief)
    push!(agb.wxp, interaction.wxp)
    push!(agb.wxm, interaction.wxm)
    push!(agb.wyp, interaction.wyp)
    push!(agb.wym, interaction.wym)
    push!(agb.dxp, interaction.dxp)
    push!(agb.dxm, interaction.dxm)
    push!(agb.dyp, interaction.dyp)
    push!(agb.dym, interaction.dym)
    return nothing
end

"""
    _beams_hits_same_shape(agb, id)

Tests if all 9 component rays at section `id` hit the same object shape.
"""
@inline function _beams_hits_same_shape(agb::AstigmaticGaussianBeamlet, id::Int)::Bool
    if id > length(agb.c.segments)
        return true
    end
    i1 = intersection(agb.c.segments[id])
    if isnothing(i1)
        for b in _aux_beams(agb)
            id > length(b.segments) && continue
            isnothing(intersection(b.segments[id])) || return false
        end
        return true
    else
        nc = normal3d(i1)
        for b in _aux_beams(agb)
            id > length(b.segments) && return false
            int = intersection(b.segments[id])
            isnothing(int) && return false
            dot(nc, normal3d(int)) > 0.5 || return false
        end
        return true
    end
end


function isparentbeam(beam::AstigmaticGaussianBeamlet, ray::AbstractRay)
    for b in _component_beams(beam)
        isparentbeam(b, ray) && return true
    end
    return false
end

"""
    _ray_to_plane_projection(plane_pos, plane_normal, ray)

Projects a ray onto a plane defined by `plane_pos` and `plane_normal`.
Returns the offset vector `h` (from plane_pos to projected point)
and the projected ray slope `u`.
"""
function _ray_to_plane_projection(plane_pos, plane_normal, ray)
    ray_pos = position(ray)
    ray_dir = direction(ray)
    d = line_plane_distance3d(plane_pos, plane_normal, ray_pos, ray_dir)
    p = ray_pos + ray_dir * d
    h = p - plane_pos
    u = ray_dir - dot(ray_dir, plane_normal) * plane_normal
    return h, u
end

"""
    parabasal_ray_parameters(agb, p0, i)

Compute the complex parabasal ray parameters (h1, u1, h2, u2) at the transverse
plane defined by the chief ray's position `p0` and direction at segment `i`.

The real parts of `h` and `u` come from the divergence rays (dxp, dyp),
while the imaginary parts come from the waist rays (wxp, wyp).
"""
function parabasal_ray_parameters(agb::AstigmaticGaussianBeamlet, p0, i)
    chief = rays(agb.c)[i]
    pn = direction(chief)

    # A Gaussian beam must model the local wavefront strictly up to its quadratic Taylor term.
    # Since ray direction is the local wavefront slope (u(x) = dW/dx = W_x + W_xx*x + 0.5*W_xxx*x^2),
    # the central difference (p - m)/2 perfectly isolates the pure paraxial curvature (W_xx)
    # and mathematically cancels out asymmetric higher-order aberrations like cubic coma (W_xxx).

    # X-axis Waist (Real Part)
    h1_real_p, u1_real_p = _ray_to_plane_projection(p0, pn, rays(agb.wxp)[i])
    h1_real_m, u1_real_m = _ray_to_plane_projection(p0, pn, rays(agb.wxm)[i])
    h1_real = (h1_real_p - h1_real_m) / 2
    u1_real = (u1_real_p - u1_real_m) / 2

    # X-axis Divergence (Imaginary Part)
    h1_imag_p, u1_imag_p = _ray_to_plane_projection(p0, pn, rays(agb.dxp)[i])
    h1_imag_m, u1_imag_m = _ray_to_plane_projection(p0, pn, rays(agb.dxm)[i])
    h1_imag = (h1_imag_p - h1_imag_m) / 2
    u1_imag = (u1_imag_p - u1_imag_m) / 2

    h1 = h1_real + im * h1_imag
    u1 = u1_real + im * u1_imag

    # Y-axis Waist (Real Part)
    h2_real_p, u2_real_p = _ray_to_plane_projection(p0, pn, rays(agb.wyp)[i])
    h2_real_m, u2_real_m = _ray_to_plane_projection(p0, pn, rays(agb.wym)[i])
    h2_real = (h2_real_p - h2_real_m) / 2
    u2_real = (u2_real_p - u2_real_m) / 2

    # Y-axis Divergence (Imaginary Part)
    h2_imag_p, u2_imag_p = _ray_to_plane_projection(p0, pn, rays(agb.dyp)[i])
    h2_imag_m, u2_imag_m = _ray_to_plane_projection(p0, pn, rays(agb.dym)[i])
    h2_imag = (h2_imag_p - h2_imag_m) / 2
    u2_imag = (u2_imag_p - u2_imag_m) / 2

    h2 = h2_real + im * h2_imag
    u2 = u2_real + im * u2_imag

    return h1, u1, h2, u2, p0
end

"""
    check_optical_invariant(agb, i; threshold = get_invariant_threshold())

Evaluate the complex optical invariant `h1 . u2 - h2 . u1 = 0` at segment `i`.
Returns `true` if the invariant holds (within `threshold`), and `false` otherwise.
"""
function check_optical_invariant(agb::AstigmaticGaussianBeamlet, i::Int;
        threshold::Real = get_invariant_threshold())
    chief = rays(agb.c)[i]
    p0 = position(chief)
    λ = wavelength(chief)
    n = refractive_index(agb, i)

    h1, u1, h2, u2, _ = parabasal_ray_parameters(agb, p0, i)

    all_pass = true

    # Check Coupling Invariant: h1 . u2 - h2 . u1 = 0
    inv_val = sum(h1 .* u2) - sum(h2 .* u1)
    if abs(inv_val) > threshold
        @warn lazy"Parabasal coupling invariant violation at segment $i: |h1.u2 - h2.u1| = $(abs(inv_val)). The paraxial astigmatic Gaussian beam tracing assumptions have broken down."
        all_pass = false
    end

    # Check Lagrange Invariant (Area Invariant): n * Im(h* . u) = λ/π
    # This ensures that each axis individually behaves like a valid Gaussian.
    H_target = λ / π
    sets = (
        (h1, u1, "x"),
        (h2, u2, "y")
    )
    for (h, u, label) in sets
        H = n * imag(sum(conj(h) .* u))
        if !isapprox(H, H_target, atol = threshold)
            @warn lazy"Lagrange invariant violation ($label) at segment $i: H=$H, target=$H_target. The beamlet has likely encountered an object that violates the paraxial assumption."
            all_pass = false
        end
    end

    return all_pass
end

"""
    parabasal_ray_parameters(agb, z)

Compute the parabasal ray parameters at distance `z` along the beam.
"""
function parabasal_ray_parameters(agb::AstigmaticGaussianBeamlet, z::Real)
    p0, i = point_on_beam(agb, z)
    return parabasal_ray_parameters(agb, p0, i)
end

"""
    shift_phase!(agb, Δϕ)

Apply a global phase shift `Δϕ` to the `AstigmaticGaussianBeamlet` by rotating the chief ray's polarization.
"""
function shift_phase!(agb::AstigmaticGaussianBeamlet, Δϕ::Real)
    for ray in rays(agb.c)
        E = polarization(ray)
        polarization!(ray, E * exp(im * Δϕ))
    end
    return nothing
end

"""
    optical_power(agb)

Compute the total integrated optical power of the beamlet. For a Gaussian beamlet,
this is typically constant through lossless propagation.
"""
function optical_power(agb::AstigmaticGaussianBeamlet)
    p0n, in_ = point_on_beam(agb, 0.0)
    chiefn = rays(agb.c)[in_]
    dirn = direction(chiefn)
    h1n, _, h2n, _, _ = parabasal_ray_parameters(agb, p0n, in_)
    area_ref = abs(_pseudo_cross2d(h1n, h2n, dirn))
    E_ref_amp = norm(polarization(chiefn))
    # Corrected power: P = I_peak * π * w0x * w0y / 2
    # We use the built-in intensity function to account for wave impedance
    I_ref = intensity(E_ref_amp)
    return 0.5 * π * area_ref * I_ref
end

"""
    waist_parameters(agb, z)

Compute the position and elliptical waist axes at distance `z` along the beam.
Returns `(p0, w1, w2, w01, w02)` where `w1` and `w2` are 3D vectors describing the semi-axes
of the beam cross-section ellipse.
"""
function waist_parameters(agb::AstigmaticGaussianBeamlet, z::Real)
    p0, i = point_on_beam(agb, z)

    # Extract the complex vectors h1 and h2 representing the amplitude matrix H from the paraxial rays.
    h1, _, h2, _, _ = parabasal_ray_parameters(agb, p0, i)

    # Physical axis estimation using the complex ray height h.
    # The true orthogonal principal axes of the beam cross-section are
    # the eigenvectors of the spot covariance matrix S.

    # Split the complex vectors into real and imaginary parts to compute the outer product.
    h1r, h1i = real(h1), imag(h1)
    h2r, h2i = real(h2), imag(h2)

    # Compute the real part of the outer tensor product H * H^dagger, yielding the spatial covariance matrix of the beam intensity.
    S = h1r * h1r' + h1i * h1i' + h2r * h2r' + h2i * h2i'

    # Perform an eigendecomposition to find the principal axes (eigenvectors) and variances (eigenvalues) of the beam cross-section.
    F = eigen(Symmetric(S))

    # F.values are sorted in ascending order. Rank is 2 (plane orthogonal to propagation direction).
    # The 3x3 matrix has rank 2 in the transverse plane, so we discard the zero eigenvalue and take the two positive variances (squared radii).
    v1 = max(zero(eltype(S)), F.values[2])
    v2 = max(zero(eltype(S)), F.values[3])

    # Scale the normalized eigenvectors by the square root of the eigenvalues to obtain the physical half-axis vectors of the elliptical spot.
    w1 = sqrt(v1) * F.vectors[:, 2]
    w2 = sqrt(v2) * F.vectors[:, 3]

    # Get waist sizes for unification
    _, _, _, _, _, w01, w02 = gauss_parameters(agb, z)

    return p0, w1, w2, w01, w02
end

"""
    _pseudo_cross2d(a, b, c)

Calculates the triple product `(a × b) ⋅ c` using a non-conjugating dot product.
"""
@inline function _pseudo_cross2d(a, b, c)
    # Scalar triple product (a × b) ⋅ c
    return (a[2] * b[3] - a[3] * b[2]) * c[1] +
           (a[3] * b[1] - a[1] * b[3]) * c[2] +
           (a[1] * b[2] - a[2] * b[1]) * c[3]
end

"""
    _pseudo_dot(a, b)

Calculates the non-conjugating dot product `a ⋅ b = Σ aᵢbᵢ`.
"""
@inline function _pseudo_dot(a, b)
    return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
end

"""
    parabasal_field(agb, r, z; E_ref_amp, area_ref, z_norm)

Compute the complex scalar electric field of the [`AstigmaticGaussianBeamlet`](@ref)
at a transverse offset `r` (a 3D vector in the plane perpendicular to the chief ray)
and longitudinal position `z`.

Uses the built-in normalization `√(area_ref / area(z))` so that the result is
physical [V/m] when `E_ref_amp` matches the initial polarization amplitude.

# Arguments

- `agb`: the astigmatic Gaussian beamlet
- `r`: transverse offset vector (must be orthogonal to the chief ray direction at `z`)
- `z`: distance along the beam
- `E_ref_amp`: reference field amplitude (auto-computed from `z_norm` if `nothing`)
- `area_ref`: reference area (auto-computed from `z_norm` if `nothing`).
- `z_norm`: longitudinal position for reference normalization (default: 0).

# Formalism
The field is computed using the **Parabasal Gaussian Beamlet** formalism:
ψ(r, z) = √(area_ref / area(z)) * exp(i * k * [z + 1/2 * rᵀ * Q(z) * r])
where area(z) = (h1 × h2) · dir is the complex beam area.
"""
function parabasal_field(
        agb::AstigmaticGaussianBeamlet,
        r::AbstractArray,
        z::Real;
        E_ref_amp::Union{Nothing, Number} = nothing,
        area_ref::Union{Nothing, Complex} = nothing,
        z_norm::Real = 0.0
)
    p0, i = point_on_beam(agb, z)
    chief = rays(agb.c)[i]
    dir = direction(chief)

    if !isorthogonal3d(dir, r; atol = Config.get_orthogonality_threshold())
        error("r must lie in plane at p0/dir (dot product: $(dot(dir, r)))")
    end

    # Lazy init of reference normalization (computed once if not supplied)
    if area_ref === nothing || E_ref_amp === nothing
        p0n, in_ = point_on_beam(agb, z_norm)
        chiefn = rays(agb.c)[in_]
        if area_ref === nothing
            dirn = direction(chiefn)
            h1n, _, h2n, _, _ = parabasal_ray_parameters(agb, p0n, in_)
            area_ref = _pseudo_cross2d(h1n, h2n, dirn)
        end
        if E_ref_amp === nothing
            # Extract complex amplitude (scalar projection) to preserve phase
            E_vec = polarization(chiefn)
            max_idx = argmax(abs.(E_vec))
            E_ref_amp = Complex(norm(E_vec) * cis(angle(E_vec[max_idx])))
        end
    end

    h1, u1, h2, u2, _ = parabasal_ray_parameters(agb, p0, i)

    area = _pseudo_cross2d(h1, h2, dir)
    # Regularization to prevent division by zero at caustics or NaN propagation
    if isnan(area) || abs(area) < 1e-25
        area = Complex(1e-25, 1e-25)
    end
    ξ1 = _pseudo_cross2d(h1, r, dir)
    ξ2 = _pseudo_cross2d(h2, r, dir)

    # Complex quadratic phase term w = r^T Q r / 2
    # We use _pseudo_dot to compute the non-conjugating dot product r^T u.
    w = (ξ1 * _pseudo_dot(u2, r) - ξ2 * _pseudo_dot(u1, r)) / (2 * area)

    # Phase correction for refractive index (OPL - geometric length)
    p_parent = agb.parent
    l_parent = isnothing(p_parent) ? 0.0 : length(p_parent)
    opl_parent = isnothing(p_parent) ? 0.0 : optical_path_length(p_parent)

    Δl = opl_parent - l_parent
    z_sum = l_parent
    for j in 1:(i - 1)
        ray_j = rays(agb.c)[j]
        Δl += optical_path_length(ray_j) - length(ray_j)
        z_sum += length(ray_j)
    end
    Δl += (refractive_index(chief) - 1) * (z - z_sum)

    k0 = 2π / wavelength(chief)

    # area_ref / area is essentially (1 / (1 + i*z/zr))^2 for stigmatic beams.
    # The sqrt gives the standard -atan(z/zr) Gouy phase natively.
    # The phase includes the OPL correction Δl to ensure coherence in media.
    ψ = sqrt(area_ref / area) * exp(im * k0 * (z + w + Δl))
    return E_ref_amp * ψ
end

"""
    electric_field(agb, r, z)

Convenience wrapper for [`parabasal_field`](@ref) using the beamlet's starting position (z=0)
as the reference normalization.
"""
function electric_field(agb::AstigmaticGaussianBeamlet, r::AbstractArray, z::Real)
    return parabasal_field(agb, r, z; z_norm = 0.0)
end

"""
    polarized_field(agb, r, z)

Compute the complex **vector** electric field [V/m] of the [`AstigmaticGaussianBeamlet`](@ref)
at position `(r, z)`. Returns a 3D vector.
"""
function polarized_field(agb::AstigmaticGaussianBeamlet, r::AbstractArray, z::Real)
    p0, i = point_on_beam(agb, z)
    chief = rays(agb.c)[i]
    E_vec = polarization(chief)
    # The complex scalar field already includes the propagation phase and Gouy phase.
    # We normalize to the chief ray's polarization magnitude to avoid double-counting amplitude.
    ψ = parabasal_field(agb, r, z; E_ref_amp = 1.0)
    return E_vec * ψ
end

"""
    intensity(agb::AstigmaticGaussianBeamlet, r, z)

Compute the optical intensity [W/m²] of the beamlet at position `(r, z)`.
"""
function intensity(agb::AstigmaticGaussianBeamlet, r::AbstractArray, z::Real)
    E = electric_field(agb, r, z)
    return intensity(E)
end

"""
    rayleigh_range(agb::AstigmaticGaussianBeamlet)

Returns the Rayleigh range for the x and y axes of the beamlet as a tuple `(z_rx, z_ry)`.
"""
function rayleigh_range(agb::AstigmaticGaussianBeamlet)
    λ = wavelength(agb)
    _, w1, w2 = waist_parameters(agb, 0.0)
    return (rayleigh_range(λ, norm(w1)), rayleigh_range(λ, norm(w2)))
end

"""
    gauss_parameters(agb::AstigmaticGaussianBeamlet, z::Real)

Compute the scalar Gaussian beam parameters at distance `z`.
Returns a tuple `(w1, w2, R1, R2, ψ, w01, w02)` where:
- `w1`, `w2`: beam radii along the principal axes
- `R1`, `R2`: radii of curvature along the principal axes
- `ψ`: total Gouy phase shift
- `w01`, `w02`: waist radii along the principal axes
"""
function gauss_parameters(agb::AstigmaticGaussianBeamlet, z::Real)
    p0, i = point_on_beam(agb, z)
    h1, u1, h2, u2, _ = parabasal_ray_parameters(agb, p0, i)
    n = refractive_index(agb, i)
    λ = wavelength(agb)
    T = typeof(λ)
    dir = direction(rays(agb.c)[i])

    # Axis 1
    w1 = norm(h1)
    # real(sum(u1 .* conj.(h1)))
    invR1 = (real(u1[1]) * real(h1[1]) + imag(u1[1]) * imag(h1[1]) +
             real(u1[2]) * real(h1[2]) + imag(u1[2]) * imag(h1[2]) +
             real(u1[3]) * real(h1[3]) + imag(u1[3]) * imag(h1[3])) / w1^2
    R1 = isapprox(invR1, 0, atol = 1e-20) ? T(Inf) : 1 / invR1
    # imag(sum(h1 .* conj.(u1)))
    H1 = abs(n * (imag(h1[1]) * real(u1[1]) - real(h1[1]) * imag(u1[1]) +
              imag(h1[2]) * real(u1[2]) - real(h1[2]) * imag(u1[2]) +
              imag(h1[3]) * real(u1[3]) - real(h1[3]) * imag(u1[3])))
    w01 = H1 / (n * norm(u1))

    # Axis 2
    w2 = norm(h2)
    # real(sum(u2 .* conj.(h2)))
    invR2 = (real(u2[1]) * real(h2[1]) + imag(u2[1]) * imag(h2[1]) +
             real(u2[2]) * real(h2[2]) + imag(u2[2]) * imag(h2[2]) +
             real(u2[3]) * real(h2[3]) + imag(u2[3]) * imag(h2[3])) / w2^2
    R2 = isapprox(invR2, 0, atol = 1e-20) ? T(Inf) : 1 / invR2
    # imag(sum(h2 .* conj.(u2)))
    H2 = abs(n * (imag(h2[1]) * real(u2[1]) - real(h2[1]) * imag(u2[1]) +
              imag(h2[2]) * real(u2[2]) - real(h2[2]) * imag(u2[2]) +
              imag(h2[3]) * real(u2[3]) - real(h2[3]) * imag(u2[3])))
    w02 = H2 / (n * norm(u2))

    # Total Gouy phase
    area = _pseudo_cross2d(h1, h2, dir)
    ψ = -0.5 * angle(area)

    return (w1, w2, R1, R2, ψ, w01, w02)
end

function waist_parameters(agb::AstigmaticGaussianBeamlet{G}, zs::AbstractArray) where {G}
    n = length(zs)
    p0 = Vector{Point3{G}}(undef, n)
    w1 = Vector{Point3{G}}(undef, n)
    w2 = Vector{Point3{G}}(undef, n)
    w01 = Vector{G}(undef, n)
    w02 = Vector{G}(undef, n)
    @inbounds for i in 1:n
        p0[i], w1[i], w2[i], w01[i], w02[i] = waist_parameters(agb, zs[i])
    end
    return (p0, w1, w2, w01, w02)
end
