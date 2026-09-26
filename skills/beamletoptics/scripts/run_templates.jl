# Smoke-test the headless templates (01-06) against the active BeamletOptics environment.
#
#   julia --project=<env with BeamletOptics> scripts/run_templates.jl
#
# 07_render_system.jl needs a Makie backend and is skipped unless `--render` is passed.
const TEMPLATE_DIR = normpath(joinpath(@__DIR__, "..", "templates"))

files = sort(filter(f -> endswith(f, ".jl"), readdir(TEMPLATE_DIR)))
"--render" in ARGS || filter!(f -> !occursin("render", f), files)

failed = String[]
for f in files
    println("\n=== ", f)
    try
        # each template runs in its own module to avoid name clashes
        m = Module(Symbol(splitext(f)[1]))
        Base.include(m, joinpath(TEMPLATE_DIR, f))
    catch err
        push!(failed, f)
        showerror(stdout, err)
        println()
    end
end

println("\n", length(files) - length(failed), "/", length(files), " templates ran without error")
isempty(failed) || (println("failed: ", join(failed, ", ")); exit(1))
