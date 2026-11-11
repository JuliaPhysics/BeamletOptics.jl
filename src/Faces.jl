"""
    isonface(face::AbstractFace, part, hit::Intersection, point, normal)

Returns `true` if the `hit` of the body `part` at the `point` with its outward `normal`, both in the
local frame of the `part`, lies on the `face`. By default, the `face` is compared with the native face
of the hit, see [`face`](@ref); a face that is defined by geometry implements its own method.
"""
isonface(face::AbstractFace, ::AbstractObject, hit::Intersection, point, normal) = face == hit.face

"""
    AllFaces()

Every face of a body, e.g. for a coating on all its surfaces.
"""
struct AllFaces <: AbstractFace end

isonface(::AllFaces, ::AbstractObject, ::Intersection, point, normal) = true

"""
    FrontFace()

The front surface of a lens: every face whose outward normal points against the optical axis of the
lens, its local y-axis, including the flat annulus of a mechanical flange. See also [`BackFace`](@ref)
and [`EdgeFace`](@ref).
"""
struct FrontFace <: AbstractFace end

"""
    BackFace()

The back surface of a lens: every face whose outward normal points along the optical axis of the lens,
its local y-axis. See also [`FrontFace`](@ref).
"""
struct BackFace <: AbstractFace end

"""
    EdgeFace()

The edge of a lens: every face whose outward normal is perpendicular to the optical axis of the lens,
its local y-axis, e.g. the cylindrical rim. See also [`FrontFace`](@ref).
"""
struct EdgeFace <: AbstractFace end

# Normals within this cosine of the transverse plane belong to the edge
const edge_cosine = 1e-6

isonface(::FrontFace, ::AbstractObject, ::Intersection, point, normal) = normal[2] < -edge_cosine
isonface(::BackFace, ::AbstractObject, ::Intersection, point, normal) = normal[2] > edge_cosine
isonface(::EdgeFace, ::AbstractObject, ::Intersection, point, normal) = abs(normal[2]) ≤ edge_cosine

"""
    FlatFace(normal)

The flat face of a body whose outward normal is `normal` in the local frame of the body, e.g. the
hypotenuse of a right-angle prism or a face of a plate.
"""
struct FlatFace{T <: Real} <: AbstractFace
    normal::Point3{T}
end

FlatFace(normal::AbstractVector{<:Real}) = FlatFace(Point3(normalize(float.(normal))...))

isonface(face::FlatFace, ::AbstractObject, ::Intersection, point, normal) =
    isapprox(normal, face.normal; atol = edge_cosine)

"""
    MaskedFace(face, mask::AbstractShape)

The part of the `face` that lies within the `mask`, a shape in the local frame of the body, e.g. a
`BeamletOptics.CylinderSDF` along the optical axis for a circular area of a lens surface, a
`BeamletOptics.RingSDF` for an annulus or a `BeamletOptics.BoxSDF` for a rectangle. The mask must
implement [`isinside`](@ref): a shape without a method contains no point, and the masked face never
applies.

```julia
center = MaskedFace(FrontFace(), BeamletOptics.CylinderSDF(5e-3, 1e-2))   # r < 5 mm
```
"""
struct MaskedFace{F <: AbstractFace, S <: AbstractShape} <: AbstractFace
    face::F
    mask::S
end

isonface(face::MaskedFace, part::AbstractObject, hit::Intersection, point, normal) =
    isonface(face.face, part, hit, point, normal) && isinside(face.mask, point)

"""
    FaceModels

The surface models of the faces of a body, `face => model` pairs in the order in which they were
set, see [`surface_model!`](@ref). The faces and models of a body are of different types.
"""
const FaceModels = Vector{Pair{AbstractFace, AbstractSurfaceModel}}

"""
    face_models(part)

Returns the surface models of the faces of the body `part`, as a vector of `face => model` pairs, see
[`surface_model!`](@ref); an empty tuple for a part whose faces carry no models.
"""
face_models(::AbstractObject) = ()

"""
    surface_model!(part, face::AbstractFace, model)

Coats the `face` of the body `part`, e.g. a lens or a prism, with the surface `model`, e.g. a
[`ThinFilmCoating`](@ref), or removes its model with `nothing`. The model acts at every boundary of
the face, between the media on both sides; the layers of a coating lie on the outside of the part.
Where several faces of the part apply to a hit, e.g. a [`MaskedFace`](@ref) and the face around it,
the face that was coated first acts.

!!! info "Cemented faces"
    Where two coated faces touch, e.g. the cemented faces of a doublet, both coatings act, each from
    its own body outward, see [`CementedCoatings`](@ref): two thin elements are cascaded, two
    [`ThinFilmCoating`](@ref)s form one stack. Other pairs, e.g. a thin-film design on one face and an
    ideal mirror on the other, cannot be combined: the coating of the body that the ray leaves acts,
    and a warning is issued once.

```julia
lens = SphericalLens(50e-3, -50e-3, 5e-3, 25e-3, ConstantRefractiveIndex(1.5))
surface_model!(lens, FrontFace(), ThinFilmCoating(QuarterWaveLayer(ConstantRefractiveIndex(1.38), 550e-9)))
```
"""
function surface_model!(part::AbstractObject, face::AbstractFace, model::Nullable{AbstractSurfaceModel})
    face_model!(face_models(part), part, face, model)
    return nothing
end

face_model!(::Tuple{}, part::AbstractObject, ::AbstractFace, model) =
    throw(ArgumentError(lazy"the faces of $(nameof(typeof(part))) carry no surface models"))

function face_model!(table::Vector, ::AbstractObject, face::AbstractFace, model::AbstractSurfaceModel)
    i = findfirst(entry -> first(entry) == face, table)
    isnothing(i) ? push!(table, face => model) : (table[i] = face => model)
    return nothing
end

function face_model!(table::Vector, ::AbstractObject, face::AbstractFace, ::Nothing)
    filter!(entry -> first(entry) != face, table)
    return nothing
end

"""
    surface_model(part, face::AbstractFace)

Returns the surface model of the `face` of the body `part`, or `nothing`, see [`surface_model!`](@ref).
"""
function surface_model(part::AbstractObject, face::AbstractFace)
    for (f, model) in face_models(part)
        f == face && return model
    end
    return nothing
end

"""
    face_model(part, hit::Intersection, point)

Returns the surface model of the face of the body `part` on which its `hit` at the `point` lies, or
`nothing`.
"""
function face_model(part::AbstractObject, hit::Intersection, point)
    table = face_models(part)
    isempty(table) && return nothing
    # Hit in the local frame of the part
    R = orientation(part)
    p = transpose(R) * (point - position(part))
    n = transpose(R) * normal3d(hit)
    for (face, model) in table
        isonface(face, part, hit, p, n) && return model
    end
    return nothing
end
