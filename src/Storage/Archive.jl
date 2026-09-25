#=
Container of a stored setup: a zip archive with a `setup.toml` description and an
`assets/` folder for binary or foreign files (meshes, large arrays, later e.g. STEP files).
=#

const SETUP_FORMAT = "BeamletOptics.setup"
const SETUP_ENTRY = "setup.toml"
const ASSET_DIR = "assets"

"Arrays with more elements than this are written as `.bin` assets instead of inline TOML arrays"
const INLINE_ARRAY_LIMIT = 1024

"""
    StorageContext

State that is shared by all [`encode`](@ref)/[`decode`](@ref) calls of one save or load.

# Fields

- `file_version`: BMO version that wrote the file (the running version when saving)
- `assets`: relative path => bytes of all assets, pending ones when saving, read on demand when loading
- `hashes`: SHA-256 hex digest => relative path, deduplicates identical assets when saving
- `reader`: the opened archive when loading, `nothing` when saving
"""
mutable struct StorageContext
    file_version::VersionNumber
    assets::Dict{String, Vector{UInt8}}
    hashes::Dict{String, String}
    reader::Union{Nothing, ZipReader}
end

StorageContext(version::VersionNumber) = StorageContext(version, Dict{String, Vector{UInt8}}(), Dict{String, String}(), nothing)

"""
    write_asset!(ctx::StorageContext, stem, ext, bytes) -> String

Adds `bytes` as the asset `assets/<stem>-<sha8>.<ext>` and returns this relative path, which is
stored in the TOML description. Identical bytes are stored once, the path of the first asset
is returned for all further ones.
"""
function write_asset!(ctx::StorageContext, stem::AbstractString, ext::AbstractString, bytes::Vector{UInt8})
    hash = bytes2hex(sha256(bytes))
    haskey(ctx.hashes, hash) && return ctx.hashes[hash]
    occursin(r"^[A-Za-z0-9_.-]+$", stem) || throw(ArgumentError("invalid asset stem \"$stem\""))
    occursin(r"^[A-Za-z0-9]+$", ext) || throw(ArgumentError("invalid asset extension \"$ext\""))
    path = "$ASSET_DIR/$stem-$(hash[1:8]).$ext"
    # Different contents with the same stem and hash prefix are practically impossible, but must not overwrite
    haskey(ctx.assets, path) && (path = "$ASSET_DIR/$stem-$hash.$ext")
    ctx.assets[path] = bytes
    ctx.hashes[hash] = path
    return path
end

"""
    read_asset(ctx::StorageContext, relpath) -> Vector{UInt8}

Returns the bytes of the asset `relpath` of the loaded archive.
"""
function read_asset(ctx::StorageContext, relpath::AbstractString)
    haskey(ctx.assets, relpath) && return ctx.assets[relpath]
    reader = ctx.reader
    if isnothing(reader) || isnothing(zip_findlast_entry(reader, relpath))
        throw(ArgumentError("asset \"$relpath\" is missing in the setup file"))
    end
    return ctx.assets[relpath] = zip_readentry(reader, relpath)
end

#=
Encoding of small values
=#

encode_vec3(v) = Float64[v[1], v[2], v[3]]
decode_vec3(x) = Point3{Float64}(x[1], x[2], x[3])

"Row-major nested vectors, i.e. `m[i][j] == M[i, j]`"
encode_mat3(M) = [Float64[M[i, 1], M[i, 2], M[i, 3]] for i in 1:3]
decode_mat3(x) = SMatrix{3, 3, Float64, 9}(x[i][j] for i in 1:3, j in 1:3)

encode_complex(z::Number) = Float64[real(z), imag(z)]
decode_complex(x) = complex(Float64(x[1]), Float64(x[2]))

const ARRAY_ELTYPES = Dict{String, DataType}(
    "Float64" => Float64, "Int64" => Int64, "ComplexF64" => ComplexF64)

"""
    encode_array(A, ctx; stem = "array") -> Dict{String, Any}

Encodes an array with element type `Float64`, `Int64` or `ComplexF64` as `eltype`, `size` and
either `data` (column-major, inline) or `asset` (little-endian `.bin` asset) if `A` has more than
$INLINE_ARRAY_LIMIT elements. Complex numbers are stored as interleaved real and imaginary parts.
"""
function encode_array(A::AbstractArray{T}, ctx::StorageContext; stem::AbstractString = "array") where {T}
    name = findfirst(==(T), ARRAY_ELTYPES)
    isnothing(name) && throw(ArgumentError("arrays of $T cannot be stored"))
    d = Dict{String, Any}("eltype" => name, "size" => collect(Int, size(A)))
    flat = collect(T, vec(A))
    if length(flat) > INLINE_ARRAY_LIMIT
        d["asset"] = write_asset!(ctx, stem, "bin", Vector{UInt8}(reinterpret(UInt8, htol.(flat))))
    else
        d["data"] = T <: Complex ? collect(reinterpret(Float64, flat)) : flat
    end
    return d
end

"""
    decode_array(d, ctx) -> Array

Inverse of [`encode_array`](@ref).
"""
function decode_array(d::AbstractDict, ctx::StorageContext)
    T = get(ARRAY_ELTYPES, d["eltype"]) do
        throw(ArgumentError("unknown array eltype \"$(d["eltype"])\""))
    end
    dims = Tuple(Int.(d["size"]))
    flat = if haskey(d, "asset")
        bytes = read_asset(ctx, d["asset"])
        length(bytes) == prod(dims) * sizeof(T) ||
            throw(ArgumentError("asset \"$(d["asset"])\" does not match size $dims of $T"))
        ltoh.(collect(reinterpret(T, bytes)))
    elseif T <: Complex
        collect(reinterpret(T, Float64.(d["data"])))
    else
        T.(d["data"])
    end
    length(flat) == prod(dims) || throw(ArgumentError("array data does not match size $dims"))
    return reshape(flat, dims)
end

#=
Archive IO
=#

# Keys that come first in each TOML table, all others follow alphabetically
const _KEY_ORDER = Dict("format" => 1, "bmo_version" => 2, "eltype" => 3, "id" => 4, "type" => 5)
_key_order(k) = (get(_KEY_ORDER, k, 100), k)

function _write_archive(path::AbstractString, setup::AbstractDict, ctx::StorageContext)
    toml = try
        sprint(io -> TOML.print(io, setup; sorted = true, by = _key_order))
    catch e
        throw(ArgumentError("setup contains values that TOML cannot represent: $(sprint(showerror, e))"))
    end
    # Write to a temporary file first, so that a failed save does not destroy an existing file
    tmp = tempname(dirname(abspath(path)); cleanup = false)
    try
        ZipWriter(tmp) do w
            zip_writefile(w, SETUP_ENTRY, codeunits(toml))
            for relpath in sort!(collect(keys(ctx.assets)))
                zip_writefile(w, relpath, ctx.assets[relpath])
            end
        end
        mv(tmp, path; force = true)
    finally
        isfile(tmp) && rm(tmp)
    end
    return path
end

function _read_archive(path::AbstractString)
    isfile(path) || throw(ArgumentError("setup file \"$path\" does not exist"))
    reader = try
        ZipReader(read(path))
    catch e
        throw(ArgumentError("\"$path\" is not a BeamletOptics setup file (no zip archive)"))
    end
    isnothing(zip_findlast_entry(reader, SETUP_ENTRY)) &&
        throw(ArgumentError("\"$path\" is not a BeamletOptics setup file ($SETUP_ENTRY is missing)"))
    setup = TOML.parse(zip_readentry(reader, SETUP_ENTRY, String))
    return setup, reader
end
