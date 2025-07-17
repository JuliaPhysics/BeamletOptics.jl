"""Ideal polarizing thin beamsplitter"""
struct PolarizingThinBeamsplitter{T<:Real,S<:AbstractShape{T}} <: AbstractBeamsplitter{T,S}
    shape::S
end

PolarizingThinBeamsplitter(width::Real, height::Real) =
    PolarizingThinBeamsplitter(RectangularFlatMesh(width, height))
PolarizingThinBeamsplitter(width::Real) = PolarizingThinBeamsplitter(width, width)

shape(bs::PolarizingThinBeamsplitter) = bs.shape

function interact3d(::AbstractSystem, bs::PolarizingThinBeamsplitter,
        beam::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray)*direction(ray)
    in_dir = direction(ray)
    s = normalize(cross(in_dir, normal))
    p = normalize(cross(s, in_dir))
    E = polarization(ray)
    Es = dot(E, s)
    Ep = dot(E, p)
    tE = Ep*p
    rE = -Es*s
    transmitted = Beam(PolarizedRay(pos, in_dir, wavelength(ray), tE))
    reflected = Beam(PolarizedRay(pos, reflection3d(in_dir, normal), wavelength(ray), rE))
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

function PolarizingCubeBeamsplitter(leg_length::Real, n::RefractiveIndex)
    front = RightAnglePrism(leg_length, leg_length, n)
    back = RightAnglePrism(leg_length, leg_length, n)
    bs = PolarizingThinBeamsplitter(√2*leg_length, leg_length)
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
