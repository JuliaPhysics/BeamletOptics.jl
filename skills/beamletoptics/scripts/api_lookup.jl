# Print docstrings and method signatures for BeamletOptics names.
#
#   julia --project=<env> scripts/api_lookup.jl SphericalLens intensity
#   julia --project=<env> scripts/api_lookup.jl --exports          # list all exported names
#   julia --project=<env> scripts/api_lookup.jl --subtypes AbstractObject
#
# Use this instead of guessing constructor arguments or keyword names.
using BeamletOptics
using InteractiveUtils: subtypes

const BMO = BeamletOptics

function docstrings(sym::Symbol)
    meta = Base.Docs.meta(BMO)
    b = Base.Docs.Binding(BMO, sym)
    haskey(meta, b) || return "(no docstring)"
    return join((join(d.text) for d in values(meta[b].docs)), "\n" * "-"^40 * "\n")
end

function lookup(name::AbstractString)
    sym = Symbol(name)
    if !isdefined(BMO, sym)
        println("`$name` is not defined in BeamletOptics.")
        hits = filter(n -> !startswith(string(n), '#') &&
                           occursin(lowercase(name), lowercase(string(n))), names(BMO; all = true))
        isempty(hits) || println("Similar names: ", join(string.(hits[1:min(end, 15)]), ", "))
        return
    end
    obj = getfield(BMO, sym)
    println("="^80)
    println(name, sym in names(BMO) ? "  (exported)" : "  (NOT exported: use BeamletOptics.$name)")
    println("="^80)
    println(docstrings(sym))
    if obj isa Function || obj isa Type
        println("-- methods --")
        for m in methods(obj)
            parentmodule(m) === BMO && println("  ", m)
        end
    end
end

if isempty(ARGS)
    println("usage: api_lookup.jl NAME [NAME...] | --exports | --subtypes TYPE")
elseif ARGS[1] == "--exports"
    foreach(println, sort(string.(names(BMO))))
elseif ARGS[1] == "--subtypes"
    walk(T, d = 0) = (println("  "^d, T); foreach(S -> walk(S, d + 1), subtypes(T)))
    walk(getfield(BMO, Symbol(ARGS[2])))
else
    foreach(lookup, ARGS)
end
