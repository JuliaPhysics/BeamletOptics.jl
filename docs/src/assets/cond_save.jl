using GLMakie

const GLOBAL_USE_PLACEHOLDERS = true

"""
    conditional_include(fname; use_placeholder=!haskey(ENV, "CI"))

This is a function to speed up compilation time when building the docs locally.
It checks whether the `include` call is made locally or from within the **CI environment**.
In the former case, script evaluation is avoided and a placeholder figure is created for 
each `save` call within the script.

Evaluation can be forced by setting `use_placeholder = false`. Global evaluation can be forced
by setting `GLOBAL_USE_PLACEHOLDERS = false` in this file (does not apply for CI env.).
"""
function conditional_include(fname::String; use_placeholder::Bool=!haskey(ENV, "CI"))
    # Check for global flag if not within CI env.
    if !haskey(ENV, "CI") && !GLOBAL_USE_PLACEHOLDERS
        use_placeholder = false
    end
    # Check if to run script
    if use_placeholder        
        content = read(fname, String)
        # match save("filename")
        pattern = r"""save\(\s*"([^"]+)"""
        filenames = [m.captures[1] for m in eachmatch(pattern, content)]
        for f in filenames
            fig = Figure(size=(600,300))
            Box(fig[1, 1], color = :gray)
            save(f, fig)
            @info "Saving placeholder for $f"
        end
    else
        @info "Running script for $fname"
        include(fname)
    end
end