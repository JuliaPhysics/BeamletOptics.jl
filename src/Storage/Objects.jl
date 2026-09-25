#=
Storage of the optical components, see `to_storage`.

Objects are stored as their type plus fields: shapes (at world pose) and materials are delegated to
`encode`/`encode_material`. Objects that consist of several sub-objects (MultiShape) store these inline.
Transient state, e.g. detector hits, is not stored.
=#

# Encodes the refractive index of `obj`, errors name the object type
function _encode_material_of(obj, n, ctx::StorageContext)
    try
        return encode_material(n, ctx)
    catch e
        e isa ArgumentError || rethrow()
        throw(ArgumentError("$(nameof(typeof(obj))): $(e.msg)"))
    end
end

#=
Objects that are defined by their shape only
=#

for (T, field) in ((Mirror, :shape), (RoundPlanoMirror, :shape), (SphericalMirror, :shape),
                   (RightAnglePrismMirror, :shape), (Retroreflector, :mesh),
                   (IntersectableObject, :shape), (NonInteractableObject, :shape))
    key = string(field)
    @eval to_storage(o::$T, ctx) = Dict{String, Any}($key => encode(getfield(o, $(QuoteNode(field))), ctx))
    @eval from_storage(::Type{$T}, d, ctx) = $T(decode(d[$key], ctx))
end

#=
Refractive optics
=#

for T in (Lens, Prism)
    @eval to_storage(o::$T, ctx) = Dict{String, Any}(
        "shape" => encode(o.shape, ctx), "n" => _encode_material_of(o, o.n, ctx))
    @eval from_storage(::Type{$T}, d, ctx) = $T(decode(d["shape"], ctx), decode_material(d["n"], ctx))
end

#=
Beamsplitters, polarizers and detectors
=#

# `reflectance` and `transmittance` are the field amplitudes as held by the struct, i.e. R² + T² = 1
to_storage(bs::ThinBeamsplitter, ctx) = Dict{String, Any}("shape" => encode(bs.shape, ctx),
    "reflectance" => Float64(bs.reflectance), "transmittance" => Float64(bs.transmittance))
from_storage(::Type{ThinBeamsplitter}, d, ctx) =
    ThinBeamsplitter(decode(d["shape"], ctx), Float64(d["reflectance"]), Float64(d["transmittance"]))

# The Jones matrix is stored as complex entries, it is real for Float64 filters
to_storage(pf::PolarizationFilter, ctx) = Dict{String, Any}("shape" => encode(pf.shape, ctx),
    "jones" => encode_array(Matrix{ComplexF64}(pf.JMat.data), ctx; stem = "jones"),
    "cutoff" => Float64(pf.cutoff))
function from_storage(::Type{PolarizationFilter}, d, ctx)
    J = decode_array(d["jones"], ctx)
    size(J) == (3, 3) || throw(ArgumentError("Jones matrix must be 3×3, got size $(size(J))"))
    all(iszero ∘ imag, J) ||
        throw(ArgumentError("complex Jones matrices are not supported by PolarizationFilter{Float64}"))
    return PolarizationFilter(decode(d["shape"], ctx), GlobalJonesBasis(SMatrix{3, 3, Float64, 9}(real.(J))),
        Float64(d["cutoff"]))
end

to_storage(det::Detector, ctx) = Dict{String, Any}("shape" => encode(det.shape, ctx), "stop" => det.stop)
from_storage(::Type{Detector}, d, ctx) = Detector(decode(d["shape"], ctx), nothing, Bool(d["stop"]), ReentrantLock())

#=
Objects made of several sub-objects, stored inline
=#

for (T, fields) in ((DoubletLens, (:front, :back)), (TripletLens, (:front, :middle, :back)),
                    (CubeBeamsplitter, (:front, :back, :coating)),
                    (RectangularPlateBeamsplitter, (:substrate, :coating)),
                    (RoundPlateBeamsplitter, (:substrate, :coating)),
                    (LinearPolarizer, (:filter, :front, :back)))
    fieldkeys = map(string, fields)
    @eval to_storage(o::$T, ctx) =
        Dict{String, Any}(k => encode(getfield(o, Symbol(k)), ctx) for k in $fieldkeys)
    @eval from_storage(::Type{$T}, d, ctx) = $T((decode(d[k], ctx) for k in $fieldkeys)...)
end

for T in (Mirror, RoundPlanoMirror, SphericalMirror, RightAnglePrismMirror, Retroreflector,
          IntersectableObject, NonInteractableObject, Lens, Prism, ThinBeamsplitter,
          PolarizationFilter, Detector, DoubletLens, TripletLens, CubeBeamsplitter,
          RectangularPlateBeamsplitter, RoundPlateBeamsplitter, LinearPolarizer)
    register_storage_type!(T, string(nameof(T)))
end
