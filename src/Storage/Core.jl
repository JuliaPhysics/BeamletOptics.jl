#=
Storage of setups, i.e. systems, sources and their pairing, see `save_setup` and `load_setup`.

Every storable type registers a tag with `register_storage_type!` and implements `to_storage` and
`from_storage`. `encode` and `decode` add or read the `"type"` key and dispatch via the registry.
=#

const STORAGE_TYPES = Dict{String, Type}()
# Keyed by the qualified type name: hashes of types are not stable across precompilation
const STORAGE_TAGS = Dict{String, String}()
_type_key(T::Type) = (w = Base.typename(T).wrapper; string(parentmodule(w), ".", nameof(w)))

"""
    register_storage_type!(T::Type, tag::String)

Registers the type `T` under the `tag` that identifies it in setup files. `T` is the type without
parameters, e.g. `Lens`, all its instances are stored under the same tag. Packages that define
storable types must call this in their `__init__` function.
"""
function register_storage_type!(T::Type, tag::String)
    wrapper = Base.typename(T).wrapper
    existing = get(STORAGE_TYPES, tag, nothing)
    if !isnothing(existing) && existing !== wrapper
        throw(ArgumentError("storage tag \"$tag\" is already registered for $existing"))
    end
    STORAGE_TYPES[tag] = wrapper
    STORAGE_TAGS[_type_key(wrapper)] = tag
    return nothing
end

"""
    storage_tag(::Type{T}) -> String

Returns the tag under which `T` is stored, see [`register_storage_type!`](@ref).
"""
function storage_tag(::Type{T}) where {T}
    get(STORAGE_TAGS, _type_key(T)) do
        throw(ArgumentError("$T cannot be stored: implement `to_storage` and `from_storage` and call `register_storage_type!`"))
    end
end

"""
    to_storage(x, ctx::StorageContext) -> Dict{String, Any}

Returns the TOML-compatible description of `x`, without the `"type"` key. Nested storable values
are described with [`encode`](@ref), large arrays with [`encode_array`](@ref).
"""
to_storage(x, ctx) = throw(ArgumentError("$(typeof(x)) cannot be stored: `to_storage` is not implemented"))

"""
    from_storage(::Type{T}, d::Dict{String, Any}, ctx::StorageContext) -> T

Rebuilds an instance of the registered type `T` from its description `d`, see [`to_storage`](@ref).
"""
from_storage(T, d, ctx) =
    throw(ArgumentError("$T cannot be loaded: `from_storage` is not implemented"))

"""
    encode(x, ctx::StorageContext) -> Dict{String, Any}

Description of `x` including its `"type"` tag.
"""
function encode(x, ctx::StorageContext)
    tag = storage_tag(typeof(x))
    d = to_storage(x, ctx)
    d["type"] = tag
    return d
end

"""
    decode(d, ctx::StorageContext)

Inverse of [`encode`](@ref).
"""
function decode(d::AbstractDict, ctx::StorageContext)
    tag = get(d, "type", nothing)
    isnothing(tag) && throw(ArgumentError("entry without \"type\": $d"))
    T = get(STORAGE_TYPES, tag) do
        throw(ArgumentError("unknown type \"$tag\" in setup written by BeamletOptics $(ctx.file_version)"))
    end
    return from_storage(T, d, ctx)
end

#=
Versioning
=#

bmo_version() = pkgversion(@__MODULE__)

_release(v::VersionNumber) = VersionNumber(v.major, v.minor, v.patch)
_series(v::VersionNumber) = v.major == 0 ? (0, Int(v.minor)) : (Int(v.major),)
_series_string(v::VersionNumber) = v.major == 0 ? "0.$(v.minor)" : "$(v.major)"

"""
    check_version(file_version; loader = bmo_version())

Throws an error if a setup written by BeamletOptics `file_version` cannot be loaded by `loader`.
Files are compatible if both versions belong to the same breaking series (same major version, or
same minor version for `0.x`) and the file is not newer than the loader.
"""
function check_version(file_version::VersionNumber; loader::VersionNumber = bmo_version())
    W, L = _release(file_version), _release(loader)
    if _series(W) != _series(L)
        throw(ArgumentError("setup written by BeamletOptics $W is incompatible with BeamletOptics $L; " *
            "load it with a BeamletOptics $(_series_string(W)).x release, e.g. `] add BeamletOptics@$(_series_string(W))`"))
    end
    W > L && throw(ArgumentError("setup written by BeamletOptics $W; update BeamletOptics to ≥ $W"))
    return nothing
end

#=
Setups
=#

# Assigns ids to systems, objects and sources
struct _IdTable
    ids::IdDict{Any, String}
    taken::Set{String}
    counts::Dict{String, Int}
end

function _IdTable(names::AbstractDict)
    table = _IdTable(IdDict{Any, String}(), Set{String}(), Dict{String, Int}())
    for (obj, name) in names
        name isa AbstractString || throw(ArgumentError("names must be strings, got $(repr(name))"))
        name in table.taken && throw(ArgumentError("duplicate name \"$name\""))
        table.ids[obj] = name
        push!(table.taken, name)
    end
    return table
end

function _id!(table::_IdTable, x, prefix::String)
    haskey(table.ids, x) && return table.ids[x]
    n = get(table.counts, prefix, 0)
    id = ""
    while true
        n += 1
        id = "$prefix-$n"
        id in table.taken || break
    end
    table.counts[prefix] = n
    push!(table.taken, id)
    return table.ids[x] = id
end

_check_eltype(x::AbstractObject{T}) where {T} =
    T === Float64 || throw(ArgumentError("only Float64 objects can be stored, got $(typeof(x))"))
_check_eltype(x) = nothing

_as_system(sys::System) = sys
_as_system(sys::StaticSystem) = System(collect(AbstractObject, sys.objects))
_as_system(sys::AbstractSystem) = throw(ArgumentError("systems of type $(typeof(sys)) cannot be stored"))

# Encodes a top-level object or source with its id, errors name the id
function _encode_entry(x, table::_IdTable, ctx)
    id = _id!(table, x, storage_tag(typeof(x)))
    entry = try
        encode(x, ctx)
    catch e
        e isa ArgumentError ? throw(ArgumentError("\"$id\": " * e.msg)) : rethrow()
    end
    entry["id"] = id
    return entry
end

# Adds `obj` and, for groups, its children to `entries` (children first), returns the id of `obj`
function _collect_object!(entries, done::IdDict{Any, String}, table::_IdTable, obj, ctx)
    haskey(done, obj) && return done[obj]
    _check_eltype(obj)
    entry = if obj isa ObjectGroup
        children = [_collect_object!(entries, done, table, child, ctx) for child in obj.objects]
        Dict{String, Any}("type" => "ObjectGroup", "id" => _id!(table, obj, "ObjectGroup"),
            "children" => children, "center" => encode_vec3(obj.center), "dir" => encode_mat3(obj.dir))
    else
        _encode_entry(obj, table, ctx)
    end
    id = entry["id"]
    push!(entries, entry)
    return done[obj] = id
end

# Sources of the right-hand side of a pair: a single source, or a tuple or vector of sources
function _pair_sources(x::Union{Tuple, AbstractVector})
    isempty(x) && throw(ArgumentError("a pair needs at least one source"))
    return collect(Any, x)
end
_pair_sources(x) = Any[x]

"""
    save_setup(path, pairs::Pair{<:AbstractSystem}...; systems = (), sources = (), names = IdDict(), metadata = Dict())

Saves the systems and sources of `pairs` to the setup file `path` (by convention `*.bmo`). A pair
is `system => source` as for `live_view`, or `system => (source1, source2, ...)` for a system that
is traced by several sources. Objects and sources that occur several times, e.g. a source that
illuminates two systems, are stored once and are shared again after [`load_setup`](@ref).

The setup file is a zip archive with a `setup.toml` description and an `assets` folder for
binary data, e.g. meshes. It records the BeamletOptics version that wrote it, see
[`load_setup`](@ref) for the compatibility rules. Only the untraced state is stored, i.e. the
sources as they were created and no detector hits.

# Arguments

- `systems`: additional systems without a source
- `sources`: additional sources without a system
- `names`: `x => "name"` for systems, objects (also inside groups) and sources, used as their ids
  in the file and returned by [`load_setup`](@ref). Other entries get ids like `"Lens-3"`.
- `metadata`: free-form dictionary with string keys and TOML-compatible values, e.g. GUI settings,
  which is stored and returned unchanged
"""
function save_setup(path::AbstractString, pairs::Pair{<:AbstractSystem}...;
        systems = (), sources = (), names::AbstractDict = IdDict{Any, String}(),
        metadata::AbstractDict = Dict{String, Any}())
    ctx = StorageContext(bmo_version())
    table = _IdTable(names)
    done = IdDict{Any, String}()
    object_entries = Dict{String, Any}[]
    system_entries = Dict{String, Any}[]
    source_entries = Dict{String, Any}[]
    add_system!(sys) = get!(done, sys) do
        stored = _as_system(sys)
        top = [_collect_object!(object_entries, done, table, obj, ctx) for obj in stored.objects]
        id = _id!(table, sys, "System")
        push!(system_entries, Dict{String, Any}("id" => id, "objects" => top))
        id
    end
    add_source!(src) = get!(done, src) do
        entry = _encode_entry(src, table, ctx)
        push!(source_entries, entry)
        entry["id"]
    end
    pair_entries = [Dict{String, Any}("system" => add_system!(first(p)),
                        "sources" => [add_source!(src) for src in _pair_sources(last(p))])
                    for p in pairs]
    foreach(add_system!, systems)
    foreach(add_source!, sources)
    setup = Dict{String, Any}(
        "format" => SETUP_FORMAT,
        "bmo_version" => string(ctx.file_version),
        "eltype" => "Float64",
        "metadata" => Dict{String, Any}(string(k) => v for (k, v) in metadata),
        "systems" => system_entries,
        "objects" => object_entries,
        "sources" => source_entries,
        "pairs" => pair_entries)
    return _write_archive(path, setup, ctx)
end

# Decodes the object `id` and, for groups, its children, detects cycles
function _resolve_object!(objects::Dict{String, Any}, entries::Dict{String, Any}, visiting::Set{String}, id, ctx)
    haskey(objects, id) && return objects[id]
    entry = get(entries, id) do
        throw(ArgumentError("unknown object id \"$id\""))
    end
    id in visiting && throw(ArgumentError("object \"$id\" contains itself"))
    push!(visiting, id)
    obj = if entry["type"] == "ObjectGroup"
        children = [_resolve_object!(objects, entries, visiting, c, ctx) for c in entry["children"]]
        group = ObjectGroup(children)
        group.center = decode_vec3(entry["center"])
        group.dir = decode_mat3(entry["dir"])
        group
    else
        decode(entry, ctx)
    end
    delete!(visiting, id)
    return objects[id] = obj
end

function _entries_by_id(list, kind)
    entries = Dict{String, Any}()
    for entry in list
        id = entry["id"]
        haskey(entries, id) && throw(ArgumentError("duplicate $kind id \"$id\""))
        entries[id] = entry
    end
    return entries
end

"""
    load_setup(path) -> (; pairs, systems, sources, names, metadata, version)

Loads a setup file written by [`save_setup`](@ref).

A file can be loaded if it was written by a BeamletOptics version of the same breaking series
(same major version, or same minor version for `0.x`) that is not newer than the running one.
Otherwise an error names the version that is needed.

# Returns

- `pairs`: `Vector` of `system => source`, or `system => (source1, source2, ...)` for pairs that
  were saved with several sources (a single source in a tuple or vector is returned unwrapped)
- `systems`: all systems, including those without a source
- `sources`: all sources, including those without a system
- `names`: id => system, object or source
- `metadata`: the `metadata` passed to [`save_setup`](@ref)
- `version`: the BeamletOptics version that wrote the file
"""
function load_setup(path::AbstractString; _loader::VersionNumber = bmo_version())
    setup, reader = _read_archive(path)
    get(setup, "format", nothing) == SETUP_FORMAT ||
        throw(ArgumentError("\"$path\" is not a BeamletOptics setup file (format is not \"$SETUP_FORMAT\")"))
    version = VersionNumber(setup["bmo_version"])
    check_version(version; loader = _loader)
    eltype = get(setup, "eltype", "Float64")
    eltype == "Float64" || throw(ArgumentError("setups with eltype $eltype are not supported"))
    ctx = StorageContext(version)
    ctx.reader = reader

    object_entries = _entries_by_id(get(setup, "objects", []), "object")
    objects = Dict{String, Any}()
    visiting = Set{String}()
    for id in keys(object_entries)
        _resolve_object!(objects, object_entries, visiting, id, ctx)
    end
    systems = Dict{String, Any}()
    system_list = System[]
    for entry in get(setup, "systems", [])
        id = entry["id"]
        (haskey(systems, id) || haskey(objects, id)) && throw(ArgumentError("duplicate id \"$id\""))
        sys = System(AbstractObject[_resolve_object!(objects, object_entries, visiting, o, ctx) for o in entry["objects"]])
        systems[id] = sys
        push!(system_list, sys)
    end
    sources = Dict{String, Any}()
    source_list = Any[]
    for entry in get(setup, "sources", [])
        id = entry["id"]
        (haskey(sources, id) || haskey(systems, id) || haskey(objects, id)) && throw(ArgumentError("duplicate id \"$id\""))
        src = decode(entry, ctx)
        sources[id] = src
        push!(source_list, src)
    end
    pairs = Pair{AbstractSystem, Any}[]
    for entry in get(setup, "pairs", [])
        sys = get(systems, entry["system"]) do
            throw(ArgumentError("pair refers to unknown system \"$(entry["system"])\""))
        end
        ids = entry["sources"]
        isempty(ids) && throw(ArgumentError("pair of system \"$(entry["system"])\" has no source"))
        srcs = map(ids) do id
            get(sources, id) do
                throw(ArgumentError("pair refers to unknown source \"$id\""))
            end
        end
        push!(pairs, sys => (length(srcs) == 1 ? only(srcs) : Tuple(srcs)))
    end
    names = merge(Dict{String, Any}(), objects, systems, sources)
    metadata = Dict{String, Any}(get(setup, "metadata", Dict{String, Any}()))
    return (; pairs, systems = system_list, sources = source_list, names, metadata, version)
end
