"""
    PolarizingThinBeamsplitter{T,S} <: AbstractBeamsplitter{T,S}

Ideal zero-thickness splitter separating s- and p-polarized components.
Reflection and transmission coefficients for s- and p-polarization can be
specified independently via the keywords `reflectance_s` and `reflectance_p`.
"""
struct PolarizingThinBeamsplitter{T<:Real,S<:AbstractShape{T}} <: AbstractBeamsplitter{T,S}
    shape::S
    rs::T
    rp::T
    ts::T
    tp::T
end

"""Placeholder shape for [`RectangularPolarizingPlateBeamsplitter`](@ref)"""
struct RectangularPolarizingPlateBeamsplitterShape{T} <: AbstractShape{T} end

"""Placeholder shape for [`RoundPolarizingPlateBeamsplitter`](@ref)"""
struct RoundPolarizingPlateBeamsplitterShape{T} <: AbstractShape{T} end

"""
    RectangularPolarizingPlateBeamsplitter <: AbstractPlateBeamsplitter

Rectangular plate beamsplitter with a polarizing coating.  The coating is
represented by a [`PolarizingThinBeamsplitter`](@ref).
"""
struct RectangularPolarizingPlateBeamsplitter{T} <: AbstractPlateBeamsplitter{T,RectangularPolarizingPlateBeamsplitterShape{T}}
    substrate::Prism{T, BoxSDF{T}}
    coating::PolarizingThinBeamsplitter{T,Mesh{T}}
end

"""
    RectangularPolarizingPlateBeamsplitter(width, height, thickness, n; reflectance_s=1.0, reflectance_p=0.0)
"""
function RectangularPolarizingPlateBeamsplitter(
        width::Real,
        height::Real,
        thickness::Real,
        n::RefractiveIndex; reflectance_s::Real=1.0, reflectance_p::Real=0.0)
    substrate_shape = BoxSDF(width, thickness, height)
    substrate = Prism(substrate_shape, n)
    translate3d!(substrate, [0, thickness/2, 0])
    coating = PolarizingThinBeamsplitter(width, height;
        reflectance_s, reflectance_p)
    zrotate3d!(coating, π)
    return RectangularPolarizingPlateBeamsplitter(substrate, coating)
end

"""
    RoundPolarizingPlateBeamsplitter(diameter, thickness, n; reflectance_s=1.0, reflectance_p=0.0)
"""
struct RoundPolarizingPlateBeamsplitter{T} <: AbstractPlateBeamsplitter{T,RoundPolarizingPlateBeamsplitterShape{T}}
    substrate::Prism{T, PlanoSurfaceSDF{T}}
    coating::PolarizingThinBeamsplitter{T,Mesh{T}}
end

function RoundPolarizingPlateBeamsplitter(
        diameter::Real,
        thickness::Real,
        n::RefractiveIndex; reflectance_s::Real=1.0, reflectance_p::Real=0.0)
    substrate_shape = PlanoSurfaceSDF(thickness, diameter)
    substrate = Prism(substrate_shape, n)
    coating = RoundPolarizingThinBeamsplitter(diameter;
        reflectance_s, reflectance_p)
    return RoundPolarizingPlateBeamsplitter(substrate, coating)
end

reflectance_s(bs::PolarizingThinBeamsplitter) = bs.rs
reflectance_p(bs::PolarizingThinBeamsplitter) = bs.rp
transmittance_s(bs::PolarizingThinBeamsplitter) = bs.ts
transmittance_p(bs::PolarizingThinBeamsplitter) = bs.tp

function PolarizingThinBeamsplitter(shape::S;
        reflectance_s::Real=1.0,
        reflectance_p::Real=0.0) where {S<:AbstractShape}
    Rs = sqrt(reflectance_s)
    Rp = sqrt(reflectance_p)
    Ts = sqrt(1 - Rs^2)
    Tp = sqrt(1 - Rp^2)
    T = promote_type(Float64, typeof(Rs), typeof(Rp))
    return PolarizingThinBeamsplitter{T,S}(shape, T(Rs), T(Rp), T(Ts), T(Tp))
end

PolarizingThinBeamsplitter(width::Real, height::Real; kw...) =
    PolarizingThinBeamsplitter(RectangularFlatMesh(width, height); kw...)
PolarizingThinBeamsplitter(width::Real; kw...) =
    PolarizingThinBeamsplitter(width, width; kw...)
RoundPolarizingThinBeamsplitter(d::Real; kw...) =
    PolarizingThinBeamsplitter(CircularFlatMesh(d/2); kw...)

shape(bs::PolarizingThinBeamsplitter) = bs.shape

function interact3d(::AbstractSystem, bs::PolarizingThinBeamsplitter,
        beam::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray)*direction(ray)
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    Jt = @SArray [transmittance_s(bs) 0 0; 0 transmittance_p(bs) 0; 0 0 1]
    Jr = @SArray [-reflectance_s(bs) 0 0; 0 reflectance_p(bs) 0; 0 0 1]
    tE = _calculate_global_E0(in_dir, in_dir, Jt, polarization(ray))
    rE = _calculate_global_E0(in_dir, out_dir, Jr, polarization(ray))
    transmitted = Beam(PolarizedRay(pos, in_dir, wavelength(ray), tE))
    reflected   = Beam(PolarizedRay(pos, out_dir, wavelength(ray), rE))
    children!(beam, [transmitted, reflected])
    return nothing
end

function interact3d(::AbstractSystem, bs::PolarizingThinBeamsplitter,
        beam::Beam{T,R}, ray::R) where {T<:Real,R<:Ray{T}}
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray)*direction(ray)
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    transmitted = Beam(Ray(pos, in_dir, wavelength(ray)))
    reflected   = Beam(Ray(pos, out_dir, wavelength(ray)))
    children!(beam, [transmitted, reflected])
    return nothing
end

"""Placeholder type for PolarizingCubeBeamsplitter"""
struct PolarizingCubeBeamsplitter{T} <: AbstractBeamsplitter{T, CubeBeamsplitterShape{T}}
    front::Prism{T, RightAnglePrismSDF{T}}
    back::Prism{T, RightAnglePrismSDF{T}}
    coating::PolarizingThinBeamsplitter{T, Mesh{T}}
end

shape_trait_of(::PolarizingCubeBeamsplitter) = MultiShape()
shape(cbs::PolarizingCubeBeamsplitter) = (cbs.front, cbs.back, cbs.coating)
refractive_index(cbs::PolarizingCubeBeamsplitter, λ::Real) = refractive_index(cbs.front, λ)

function PolarizingCubeBeamsplitter(
        leg_length::Real,
        n::RefractiveIndex;
        reflectance_s::Real=1.0,
        reflectance_p::Real=0.0)
    front = RightAnglePrism(leg_length, leg_length, n)
    back = RightAnglePrism(leg_length, leg_length, n)
    bs = PolarizingThinBeamsplitter(√2*leg_length, leg_length;
        reflectance_s, reflectance_p)
    zrotate3d!(back, deg2rad(180))
    zrotate3d!(bs, deg2rad(180-45))
    set_new_origin3d!(shape(bs))
    return PolarizingCubeBeamsplitter(front, back, bs)
end

function interact3d(system::AbstractSystem, cbs::PolarizingCubeBeamsplitter,
        beam::Beam{T,R}, ray::R) where {T<:Real,R<:AbstractRay{T}}
    if shape(intersection(ray)) === shape(cbs.front)
        interaction = interact3d(system, cbs.front, beam, ray)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
    if shape(intersection(ray)) === shape(cbs.coating)
        interact3d(system, cbs.coating, beam, ray)
        _n = refractive_index(cbs, wavelength(ray))
        refractive_index!(first(rays(beam.children[1])), _n)
        refractive_index!(first(rays(beam.children[2])), _n)
        return nothing
    end
    if shape(intersection(ray)) === shape(cbs.back)
        interaction = interact3d(system, cbs.back, beam, ray)
        hint!(interaction, Hint(cbs, shape(cbs.coating)))
        return interaction
    end
end
