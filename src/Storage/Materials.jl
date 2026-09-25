#=
Storage of refractive indices, see `encode_material` and `register_material!`.

Parametric refractive indices (`ConstantRefractiveIndex`, `SellmeierEquation`, `DiscreteRefractiveIndex`)
are stored with their parameters. Other callables, e.g. closures, are stored by the name under which
they were registered with `register_material!` and must be registered again before loading.
=#

"name => refractive index, see [`register_material!`](@ref)"
const MATERIALS = Dict{String, Any}()

"""
    register_material!(name::AbstractString, n)

Registers the [`RefractiveIndex`](@ref) `n`, e.g. a function `λ -> n`, under `name`, so that objects
using `n` can be stored with [`save_setup`](@ref). The setup file stores only the `name`, hence the
same `name` must be registered before the file is loaded with [`load_setup`](@ref), e.g.:

```julia
my_glass(λ) = 1.5 + 1e-15 / λ^2
register_material!("my_glass", my_glass)
```

Registering a `name` again replaces its refractive index.
[`ConstantRefractiveIndex`](@ref), [`SellmeierEquation`](@ref) and [`DiscreteRefractiveIndex`](@ref)
are stored with their parameters and do not need to be registered.
"""
function register_material!(name::AbstractString, n)
    isempty(name) && throw(ArgumentError("material name must not be empty"))
    test_refractive_index_function(n)
    MATERIALS[String(name)] = n
    return nothing
end

# Name of the registered material `n`, `nothing` if it is not registered
function _material_name(n)
    for name in sort!(collect(keys(MATERIALS)))
        MATERIALS[name] === n && return name
    end
    return nothing
end

_is_storable_type(T::Type) = haskey(STORAGE_TAGS, _type_key(T))

"""
    encode_material(n, ctx::StorageContext) -> Dict{String, Any}

Describes the refractive index `n`: storable types (e.g. [`ConstantRefractiveIndex`](@ref)) with
their parameters, other callables by the name they were registered with via [`register_material!`](@ref).
Throws an `ArgumentError` if `n` is neither.
"""
function encode_material(n, ctx::StorageContext)
    _is_storable_type(typeof(n)) && return encode(n, ctx)
    name = _material_name(n)
    isnothing(name) || return Dict{String, Any}("type" => "Named", "name" => name)
    throw(ArgumentError("refractive index of type $(typeof(n)) cannot be stored: " *
        "register it with `register_material!(name, n)`, or use a ConstantRefractiveIndex, " *
        "SellmeierEquation or DiscreteRefractiveIndex"))
end

"""
    decode_material(d, ctx::StorageContext)

Inverse of [`encode_material`](@ref). Named materials must be registered with [`register_material!`](@ref).
"""
function decode_material(d::AbstractDict, ctx::StorageContext)
    get(d, "type", nothing) == "Named" || return decode(d, ctx)
    name = d["name"]
    return get(MATERIALS, name) do
        throw(ArgumentError("material \"$name\" is not registered: call `register_material!(\"$name\", n)` before loading"))
    end
end

#=
Parametric refractive indices
=#

to_storage(c::ConstantRefractiveIndex, ctx) = Dict{String, Any}("n" => Float64(c.n))
from_storage(::Type{ConstantRefractiveIndex}, d, ctx) = ConstantRefractiveIndex(Float64(d["n"]))

to_storage(s::SellmeierEquation, ctx) = Dict{String, Any}(
    "B" => Float64[s.B1, s.B2, s.B3], "C" => Float64[s.C1, s.C2, s.C3])
function from_storage(::Type{SellmeierEquation}, d, ctx)
    B, C = Float64.(d["B"]), Float64.(d["C"])
    (length(B) == 3 && length(C) == 3) ||
        throw(ArgumentError("Sellmeier material needs three B and three C coefficients"))
    return SellmeierEquation(B..., C...)
end

function to_storage(dri::DiscreteRefractiveIndex, ctx)
    λs = sort!(collect(keys(dri.data)))
    return Dict{String, Any}(
        "lambda" => encode_array(Float64.(λs), ctx; stem = "lambda"),
        "n" => encode_array(Float64[dri.data[λ] for λ in λs], ctx; stem = "n"))
end
from_storage(::Type{DiscreteRefractiveIndex}, d, ctx) =
    DiscreteRefractiveIndex(decode_array(d["lambda"], ctx), decode_array(d["n"], ctx))

register_storage_type!(ConstantRefractiveIndex, "Constant")
register_storage_type!(SellmeierEquation, "Sellmeier")
register_storage_type!(DiscreteRefractiveIndex, "Discrete")
