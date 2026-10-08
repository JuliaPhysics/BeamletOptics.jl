#=
Human-readable properties of objects, shapes and beams for display, e.g. in the inspector of
the app layout of the GUI package BeamletOpticsGUI
=#

"""
    properties(x) -> Vector{Pair{String, Any}}

Returns the properties of `x` (an [`AbstractObject`](@ref), shape, beam or source) as a list of
`name => value` pairs for display, e.g. in the inspector of the GUI package BeamletOpticsGUI. Only stored or
cheaply derived values are listed, e.g. the reflectance of a beamsplitter (the power ratio R, as
passed to its constructor) or the number of hits of a [`Detector`](@ref), no design parameters
that the object does not store. The list is computed on each call, i.e. it shows the current pose
and hits.

# Conventions

- values are numbers in SI units, vectors of them (e.g. a position), `String`s or `Bool`s
- the unit of a value is given in brackets at the end of its name, e.g. `"Diameter [m]"` or
  `"Angle [rad]"`; names without brackets are dimensionless, e.g. `"Reflectance"`
- the list starts with [`default_properties`](@ref): `"Type"`, the pose and, for objects, the
  shape

```julia
julia> properties(SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3))
8-element Vector{Pair{String, Any}}:
          "Type" => "Lens"
  "Position [m]" => [0.0, 0.0, 0.0]
  "Optical axis" => [0.0, 1.0, 0.0]
         "Shape" => "UnionSDF"
 "Thickness [m]" => 0.005
   "Index model" => "function of λ"
         "n(λ₀)" => 1.5
        "λ₀ [m]" => 1.0e-6
```

The refractive index is evaluated at λ₀ = [`get_default_wavelength`](@ref)`()`. The radii and the
diameter of a lens are not listed, since its shape does not store them.

# Custom components

Add a method for your type that extends the default list, e.g.

```julia
BeamletOptics.properties(x::MyFilter) =
    [BeamletOptics.default_properties(x); "Optical density" => x.od; "Center wavelength [m]" => x.λc]
```

Without a method, [`default_properties`](@ref) is shown.
"""
properties(x) = default_properties(x)

"""
    default_properties(x) -> Vector{Pair{String, Any}}

Returns the properties that [`properties`](@ref) lists for any `x`, see its conventions:

- `"Type"`: the name of the type of `x`
- depending on the [`BeamletOptics.AbstractKinematicTrait`](@ref) of `x`: `"Position [m]"` and
  `"Optical axis"` (local +y axis) for oriented objects, `"Position [m]"` and `"Direction"` for
  directed ones (beams), only `"Position [m]"` for static objects
- for an [`AbstractObject`](@ref), depending on its [`BeamletOptics.AbstractShapeTrait`](@ref):
  the [`properties`](@ref) of its shape (`"Shape"` and the stored dimensions, e.g.
  `"Diameter [m]"`) or the number of its `"Parts"`
- for an object group: the number of its objects as `"Parts"`

Methods of `properties` for specific types start from this list, see [`properties`](@ref).
"""
default_properties(x) = Pair{String, Any}["Type" => _type_name(x); default_properties(kinematic_trait_of(x), x)]

default_properties(x::AbstractObject) = Pair{String, Any}["Type" => _type_name(x);
    default_properties(kinematic_trait_of(x), x); default_properties(shape_trait_of(x), x)]

default_properties(x::AbstractObjectGroup) = Pair{String, Any}["Type" => _type_name(x);
    default_properties(kinematic_trait_of(x), x); "Parts" => length(objects(x))]

# Pose by the kinematic trait
default_properties(::Static, _) = Pair{String, Any}[]
default_properties(::Static, x::ObjectOrGroup) = Pair{String, Any}["Position [m]" => _vector(position(x))]
default_properties(::Movable{Oriented}, x) = Pair{String, Any}["Position [m]" => _vector(position(x)),
    "Optical axis" => _vector(direction(x))]
default_properties(::Movable{Directed}, x) = Pair{String, Any}["Position [m]" => _vector(position(x)),
    "Direction" => _vector(direction(x))]

# Shape by the shape trait
default_properties(::SingleShape, x::AbstractObject) = properties(shape(x))
default_properties(::MultiShape, x::AbstractObject) = Pair{String, Any}["Parts" => length(shape(x))]

_type_name(x) = string(nameof(typeof(x)))
_vector(p) = collect(Float64, p)

#=
Shapes: the type and the stored dimensions
=#

properties(s::AbstractShape) = Pair{String, Any}["Shape" => _type_name(s)]
properties(m::AbstractMesh) = Pair{String, Any}["Shape" => _type_name(m), "Faces" => size(faces(m), 1)]
properties(s::AbstractLensSDF) = Pair{String, Any}["Shape" => _type_name(s),
    "Diameter [m]" => diameter(s), "Thickness [m]" => thickness(s)]
properties(s::AbstractSphericalSurfaceSDF) = Pair{String, Any}["Shape" => _type_name(s),
    "Radius [m]" => radius(s), "Diameter [m]" => diameter(s), "Thickness [m]" => thickness(s)]
properties(s::Union{UnionSDF, DifferenceSDF}) = Pair{String, Any}["Shape" => _type_name(s),
    "Thickness [m]" => thickness(s)]
properties(s::BoxSDF) = Pair{String, Any}["Shape" => _type_name(s), "Size [m]" => _vector(2 .* s.dimensions)]
properties(s::ConicSDF) = Pair{String, Any}["Shape" => _type_name(s), "Diameter [m]" => s.diameter,
    "Thickness [m]" => s.thickness, "Conic constant" => s.k]

"""
    _flat_size(s::AbstractShape) -> Vector{Pair{String, Any}}

Returns `"Size [m]"`, the extent of the flat mesh `s` along its local x and z axes, i.e. the width
and height of e.g. a [`Detector`](@ref), or nothing for other shapes. Linear in the number of
vertices, hence only used for flat components.
"""
_flat_size(::AbstractShape) = Pair{String, Any}[]
function _flat_size(m::AbstractMesh)
    local_vertices = (vertices(m) .- position(m)') * orientation(m)
    widths = [-(reverse(extrema(view(local_vertices, :, k)))...) for k in (1, 3)]
    return Pair{String, Any}["Size [m]" => widths]
end

#=
Components
=#

"""
Returns the refractive index model `n` (see [`RefractiveIndex`](@ref)) and its value at the
default wavelength λ₀, see [`get_default_wavelength`](@ref).
"""
function _index_properties(n)
    λ0 = get_default_wavelength()
    return Pair{String, Any}["Index model" => _index_model(n), "n(λ₀)" => n(λ0), "λ₀ [m]" => λ0]
end
# Only defined at its wavelengths, which are listed instead of a value at λ₀
function _index_properties(n::DiscreteRefractiveIndex)
    λs = sort!(collect(keys(n.data)))
    return Pair{String, Any}["Index model" => _index_model(n), "Wavelengths [m]" => λs,
        "Indices" => [n.data[λ] for λ in λs]]
end
_index_model(::Function) = "function of λ"
_index_model(n) = _type_name(n)

properties(x::AbstractRefractiveOptic) = [default_properties(x); _index_properties(x.n)]

properties(x::Union{DoubletLens, TripletLens}) =
    [default_properties(x); "Thickness [m]" => thickness(x)]

"""
Returns the power `"Reflectance"` R and `"Transmittance"` T of the `coating`, which stores the
amplitude factors √R and √T, i.e. R as passed to the constructors of the beamsplitters.
"""
_splitting_ratio(coating::ThinBeamsplitter) = Pair{String, Any}[
    "Reflectance" => reflectance(coating)^2, "Transmittance" => transmittance(coating)^2]

properties(bs::ThinBeamsplitter) = [default_properties(bs); _flat_size(shape(bs)); _splitting_ratio(bs)]
properties(bs::AbstractPlateBeamsplitter) = [default_properties(bs); _splitting_ratio(coating(bs));
    "Thickness [m]" => thickness(substrate(bs)); _index_properties(substrate(bs).n)]
properties(bs::CubeBeamsplitter) =
    [default_properties(bs); _splitting_ratio(bs.coating); _index_properties(bs.front.n)]

properties(pf::PolarizationFilter) = [default_properties(pf); _flat_size(shape(pf));
    "Transmission axis" => _vector(transmission_axis(pf))]
properties(lp::LinearPolarizer) = [default_properties(lp);
    "Transmission axis" => _vector(transmission_axis(lp)); "Thickness [m]" => thickness(lp);
    _index_properties(lp.front.n)]

"""
    hit_count(d::Detector) -> Int

The number of hits stored in the detector `d`, `0` before the first solve.
"""
hit_count(d::Detector) = isnothing(hits(d)) ? 0 : length(hits(d))

properties(d::Detector) = [default_properties(d); _flat_size(shape(d)); "Hits" => hit_count(d);
    "Stops beams" => stop(d)]

#=
Beams and sources
=#

properties(b::Beam) = [default_properties(b); "Wavelength [m]" => wavelength(first(rays(b)))]
properties(b::GaussianBeamlet) = [default_properties(b); "Wavelength [m]" => wavelength(b);
    "Waist radius [m]" => beam_waist(b); "Rayleigh range [m]" => rayleigh_range(b)]
properties(b::AstigmaticGaussianBeamlet) = [default_properties(b); "Wavelength [m]" => wavelength(b)]
properties(bg::AbstractBeamGroup) =
    [default_properties(bg); "Beams" => length(bg); "Wavelength [m]" => wavelength(bg)]
properties(cs::CollimatedSource) = [default_properties(cs); "Beams" => length(cs);
    "Wavelength [m]" => wavelength(cs); "Diameter [m]" => diameter(cs)]
properties(ps::PointSource) = [default_properties(ps); "Beams" => length(ps);
    "Wavelength [m]" => wavelength(ps); "Numerical aperture" => numerical_aperture(ps)]
