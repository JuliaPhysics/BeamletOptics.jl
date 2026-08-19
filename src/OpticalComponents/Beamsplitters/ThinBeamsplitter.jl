"""
    ThinBeamsplitter <: AbstractBeamsplitter

Represents a 2D beam-splitting device.

# Fields

- `shape`: 2D [`AbstractShape`](@ref) at which the splitting process occurs (e.g. a 2D-[`Mesh`](@ref))
- `reflectance`: scalar reflection factor
- `transmittance`: scalar transmission factor

!!! warning
    Note that the `transmittance` should be calculated from an input `reflectance` in order
    to ensure that R² + T² = 1.
"""
struct ThinBeamsplitter{T <: Real, S <: AbstractShape{T}} <: AbstractBeamsplitter{T}
    shape::S
    reflectance::T
    transmittance::T
end

reflectance(bs::ThinBeamsplitter) = bs.reflectance
transmittance(bs::ThinBeamsplitter) = bs.transmittance

is_thin_interface(::ThinBeamsplitter) = true

"""
    ThinBeamsplitter(width, height; reflectance=0.5)

Creates a zero-thickness, lossless, non-polarizing 2D rectangular [`ThinBeamsplitter`](@ref) where

- `width`: is the x-dir. edge length in [m]
- `height`: is the z-dir. edge length in [m]
- `reflectance`: kw-arg that determines how much light is **reflected**, i.e. 0.7 for a 70:30 splitter

# Additional information

!!! info "Reflectance"
    The input value for the `reflectance` R is normed such that R² + T² = 1, where T is the `transmittance`.
    The transmittance is calculated via T = √(1 - R²).

!!! warning "Reflection phase jump"
    Note that the reflection phase jump θᵣ is implemented by the individual [`interact3d`](@ref)-methods. Refer to them for more information.
"""
function ThinBeamsplitter(width::Real, height::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = RectangularFlatMesh(width, height)
    Reflected = sqrt(reflectance)
    Transmitted = sqrt(1 - Reflected^2)
    return ThinBeamsplitter(shape, Reflected, Transmitted)
end

"""
    RoundThinBeamsplitter(diameter; reflectance=0.5)

Creates a zero-thickness, 2D round [`ThinBeamsplitter`](@ref) with the specified `diameter` in [m].
For more information, refer to the [`ThinBeamsplitter`](@ref) constructor.
"""
function RoundThinBeamsplitter(diameter::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = CircularFlatMesh(diameter / 2)
    R = sqrt(reflectance)
    T = sqrt(1 - R^2)
    return ThinBeamsplitter(shape, R, T)
end

function ThinBeamsplitter(width::Real; reflectance::Real = 0.5)
    ThinBeamsplitter(width, width; reflectance)
end

Base.isvalid(bs::AbstractBeamsplitter) = reflectance(bs)^2 + transmittance(bs)^2 ≈ 1

@inline function _beamsplitter_media(system::AbstractSystem, ray::AbstractRay, int::Intersection)
    λ = wavelength(ray)
    n_incident = refractive_index(ray)
    n_sys = refractive_index(system, λ)
    
    n_exit = isnothing(int.coincident_object) ? n_sys : refractive_index(int.coincident_object, λ)
    n_enter = isnothing(int.coincident_object_2) ? n_sys : refractive_index(int.coincident_object_2, λ)
    
    tol = get_index_matching_tolerance()
    if abs(n_incident - n_exit) < tol
        n_reflected = n_exit
        n_transmitted = n_enter
    else
        n_reflected = n_enter
        n_transmitted = n_exit
    end
    return n_transmitted, n_reflected
end

@inline function _beamsplitter_transmitted_beam(
        ::AbstractBeamsplitter, ::Beam{T, R}, ray::R, n::Real, dir::AbstractVector) where {T <: Real, R <: Ray{T}}
    pos = position(ray) + length(ray) * direction(ray)
    return Beam(Ray{T}(Point3{T}(pos), Point3{T}(dir), nothing, wavelength(ray), n))
end

@inline function _beamsplitter_reflected_beam(
        ::AbstractBeamsplitter, ::Beam{T, R}, ray::R, n::Real, dir::AbstractVector) where {T <: Real, R <: Ray{T}}
    pos = position(ray) + length(ray) * direction(ray)
    return Beam(Ray{T}(Point3{T}(pos), Point3{T}(dir), nothing, wavelength(ray), n))
end

@inline function _beamsplitter_transmitted_beam(bs::AbstractBeamsplitter, ::Beam{T, R},
        ray::R, n::Real, dir::AbstractVector) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis(transmittance(bs), 0, 0, transmittance(bs))
    pos = position(ray) + length(ray) * direction(ray)
    E0 = _calculate_global_E0(bs, ray, dir, J)
    return Beam(PolarizedRay{T}(pos, dir, nothing, wavelength(ray), n, E0))
end

@inline function _beamsplitter_reflected_beam(bs::AbstractBeamsplitter, ::Beam{T, R},
        ray::R, n::Real, dir::AbstractVector) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis(-reflectance(bs), 0, 0, reflectance(bs))
    pos = position(ray) + length(ray) * direction(ray)
    E0 = _calculate_global_E0(bs, ray, dir, J)
    return Beam(PolarizedRay{T}(pos, dir, nothing, wavelength(ray), n, E0))
end

function interact3d(system::AbstractSystem, bs::ThinBeamsplitter, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    int = intersection(ray)
    n_transmitted, n_reflected = _beamsplitter_media(system, ray, int)
    
    normal = normal3d(int)
    if dot(direction(ray), normal) > 0
        normal = -normal
    end
    dir_t, _ = refraction3d(direction(ray), normal, refractive_index(ray), n_transmitted)
    dir_r = reflection3d(direction(ray), normal)
    
    children!(beam,
        (_beamsplitter_transmitted_beam(bs, beam, ray, n_transmitted, dir_t),
         _beamsplitter_reflected_beam(bs, beam, ray, n_reflected, dir_r)))
    return nothing
end

@inline function _beamsplitter_transmitted_beam(
        bs::AbstractBeamsplitter, gauss::GaussianBeamlet, ray_id::Int, n::Real, c_dir::AbstractVector, w_dir::AbstractVector, d_dir::AbstractVector)
    chief = _beamsplitter_transmitted_beam(bs, gauss.chief, rays(gauss.chief)[ray_id], n, c_dir)
    waist = _beamsplitter_transmitted_beam(bs, gauss.waist, rays(gauss.waist)[ray_id], n, w_dir)
    divergence = _beamsplitter_transmitted_beam(bs, gauss.divergence, rays(gauss.divergence)[ray_id], n, d_dir)
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = transmittance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

@inline function _beamsplitter_reflected_beam(
        bs::AbstractBeamsplitter, gauss::GaussianBeamlet, ray_id::Int, n::Real, c_dir::AbstractVector, w_dir::AbstractVector, d_dir::AbstractVector)
    chief = _beamsplitter_reflected_beam(bs, gauss.chief, rays(gauss.chief)[ray_id], n, c_dir)
    waist = _beamsplitter_reflected_beam(bs, gauss.waist, rays(gauss.waist)[ray_id], n, w_dir)
    divergence = _beamsplitter_reflected_beam(bs, gauss.divergence, rays(gauss.divergence)[ray_id], n, d_dir)
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = reflectance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

function interact3d(
        system::AbstractSystem, bs::ThinBeamsplitter, gauss::GaussianBeamlet, ray_id::Int)
    c_ray = gauss.chief.rays[ray_id]
    w_ray = gauss.waist.rays[ray_id]
    d_ray = gauss.divergence.rays[ray_id]
    
    int = intersection(c_ray)
    n_transmitted, n_reflected = _beamsplitter_media(system, c_ray, int)
    
    normal = normal3d(int)
    if dot(direction(c_ray), normal) > 0
        normal = -normal
    end
    
    c_dir_t, _ = refraction3d(direction(c_ray), normal, refractive_index(c_ray), n_transmitted)
    w_dir_t, _ = refraction3d(direction(w_ray), normal, refractive_index(w_ray), n_transmitted)
    d_dir_t, _ = refraction3d(direction(d_ray), normal, refractive_index(d_ray), n_transmitted)
    
    c_dir_r = reflection3d(direction(c_ray), normal)
    w_dir_r = reflection3d(direction(w_ray), normal)
    d_dir_r = reflection3d(direction(d_ray), normal)
    
    t = _beamsplitter_transmitted_beam(bs, gauss, ray_id, n_transmitted, c_dir_t, w_dir_t, d_dir_t)
    r = _beamsplitter_reflected_beam(bs, gauss, ray_id, n_reflected, c_dir_r, w_dir_r, d_dir_r)
    
    df = dot(direction(c_ray), normal3d(int))
    ϕ = df < 0 ? π : 0
    electric_field!(r, electric_field(r) * exp(im * ϕ))

    children!(gauss, (t, r))
    return nothing
end

@inline function _beamsplitter_transmitted_beam(
        bs::AbstractBeamsplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int, n::Real, dirs::Tuple)
    c = _beamsplitter_transmitted_beam(bs, agb.c, rays(agb.c)[ray_id], n, dirs[1])
    wxp = _beamsplitter_transmitted_beam(bs, agb.wxp, rays(agb.wxp)[ray_id], n, dirs[2])
    wxm = _beamsplitter_transmitted_beam(bs, agb.wxm, rays(agb.wxm)[ray_id], n, dirs[3])
    wyp = _beamsplitter_transmitted_beam(bs, agb.wyp, rays(agb.wyp)[ray_id], n, dirs[4])
    wym = _beamsplitter_transmitted_beam(bs, agb.wym, rays(agb.wym)[ray_id], n, dirs[5])
    dxp = _beamsplitter_transmitted_beam(bs, agb.dxp, rays(agb.dxp)[ray_id], n, dirs[6])
    dxm = _beamsplitter_transmitted_beam(bs, agb.dxm, rays(agb.dxm)[ray_id], n, dirs[7])
    dyp = _beamsplitter_transmitted_beam(bs, agb.dyp, rays(agb.dyp)[ray_id], n, dirs[8])
    dym = _beamsplitter_transmitted_beam(bs, agb.dym, rays(agb.dym)[ray_id], n, dirs[9])
    return AstigmaticGaussianBeamlet(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym)
end

@inline function _beamsplitter_reflected_beam(
        bs::AbstractBeamsplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int, n::Real, dirs::Tuple)
    c = _beamsplitter_reflected_beam(bs, agb.c, rays(agb.c)[ray_id], n, dirs[1])
    wxp = _beamsplitter_reflected_beam(bs, agb.wxp, rays(agb.wxp)[ray_id], n, dirs[2])
    wxm = _beamsplitter_reflected_beam(bs, agb.wxm, rays(agb.wxm)[ray_id], n, dirs[3])
    wyp = _beamsplitter_reflected_beam(bs, agb.wyp, rays(agb.wyp)[ray_id], n, dirs[4])
    wym = _beamsplitter_reflected_beam(bs, agb.wym, rays(agb.wym)[ray_id], n, dirs[5])
    dxp = _beamsplitter_reflected_beam(bs, agb.dxp, rays(agb.dxp)[ray_id], n, dirs[6])
    dxm = _beamsplitter_reflected_beam(bs, agb.dxm, rays(agb.dxm)[ray_id], n, dirs[7])
    dyp = _beamsplitter_reflected_beam(bs, agb.dyp, rays(agb.dyp)[ray_id], n, dirs[8])
    dym = _beamsplitter_reflected_beam(bs, agb.dym, rays(agb.dym)[ray_id], n, dirs[9])
    return AstigmaticGaussianBeamlet(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym)
end

function interact3d(
        system::AbstractSystem, bs::ThinBeamsplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    c_ray = rays(agb.c)[ray_id]
    int = intersection(c_ray)
    n_transmitted, n_reflected = _beamsplitter_media(system, c_ray, int)
    
    normal = normal3d(int)
    if dot(direction(c_ray), normal) > 0
        normal = -normal
    end
    
    dirs_t = (
        refraction3d(rays(agb.c)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wxp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wxm)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wyp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wym)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dxp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dxm)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dyp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dym)[ray_id], n_transmitted)[1]
    )
    dirs_r = (
        reflection3d(direction(rays(agb.c)[ray_id]), normal3d(intersection(rays(agb.c)[ray_id]))),
        reflection3d(direction(rays(agb.wxp)[ray_id]), normal3d(intersection(rays(agb.wxp)[ray_id]))),
        reflection3d(direction(rays(agb.wxm)[ray_id]), normal3d(intersection(rays(agb.wxm)[ray_id]))),
        reflection3d(direction(rays(agb.wyp)[ray_id]), normal3d(intersection(rays(agb.wyp)[ray_id]))),
        reflection3d(direction(rays(agb.wym)[ray_id]), normal3d(intersection(rays(agb.wym)[ray_id]))),
        reflection3d(direction(rays(agb.dxp)[ray_id]), normal3d(intersection(rays(agb.dxp)[ray_id]))),
        reflection3d(direction(rays(agb.dxm)[ray_id]), normal3d(intersection(rays(agb.dxm)[ray_id]))),
        reflection3d(direction(rays(agb.dyp)[ray_id]), normal3d(intersection(rays(agb.dyp)[ray_id]))),
        reflection3d(direction(rays(agb.dym)[ray_id]), normal3d(intersection(rays(agb.dym)[ray_id])))
    )

    t = _beamsplitter_transmitted_beam(bs, agb, ray_id, n_transmitted, dirs_t)
    r = _beamsplitter_reflected_beam(bs, agb, ray_id, n_reflected, dirs_r)

    df = dot(direction(c_ray), normal3d(int))
    ϕ = df < 0 ? π : 0
    # Add conditional phase flip to reflected beam chief polarization
    chief_ray = first(rays(r.c))
    E_new = polarization(chief_ray) * exp(im * ϕ)
    polarization!(chief_ray, E_new)

    children!(agb, (t, r))
    return nothing
end
