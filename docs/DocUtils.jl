module DocUtils

using GLMakie: Figure, Axis, hidedecorations!, text!, save
using Dates: now

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
        # match save("filename") but ignore comments, i.e. # save(...)
        pattern = Regex("^[ \\t]*save\\(\\s*\"([^\"]+)\"\\s*,\\s*([a-zA-Z_][a-zA-Z0-9_]*)", "m")
        for match in eachmatch(pattern, content)
            save_name = match.captures[1]
            fig_name = match.captures[2]
            offset = match.offset
            prior_content = content[1:offset]
            # find last missing fig_name variable match
            npattern = Regex("\\Q$fig_name\\E\\s*=\\s*Figure\\(.*?size\\s*=\\s*\\(([\\d]+)\\s*,\\s*([\\d]+)\\)")
            matches = collect(eachmatch(npattern, prior_content))
            if !isempty(matches)
                nmatch = last(matches)
                _width =  parse(Int, nmatch.captures[1])
                _height = parse(Int, nmatch.captures[2])
            else
                @info "Using default height for $fig_name"
                _width =  600
                _height = 300
            end
            # save placeholder
            fig = Figure(size=(_width, _height))
            ax = Axis(fig[1,1], backgroundcolor=:gray)
            hidedecorations!(ax)
            text!(ax, 0, 0;
                text="$save_name\n($_width x $_height)",
                align=(:center, :center),
                color=:white,
                fontsize=20
            )
            save(save_name, fig)
            @info "Saving placeholder for $save_name (Size: $_width x $_height, fig_var_name=$fig_name)"
        end
    else
        @info "Running script for $fname"
        t1 = now()
        include(fname)
        t2 = now()
        @info " Runtime: $(t2-t1)"
    end
end

function replace_build_with_src(path::String)
    clean_path = abspath(path)
    build_segment = joinpath("docs", "build")
    src_segment   = joinpath("docs", "src")
    src_path = replace(clean_path, build_segment => src_segment)
    return src_path
end

"""
    prerender_include(fname::String, cname::String)


"""
function prerender_include(fname::String, cname::String)
    location = replace_build_with_src(dirname(cname))
    content = read(fname, String)
    pattern = Regex("^[ \\t]*save\\(\\s*\"([^\"]+)\"\\s*,\\s*([a-zA-Z_][a-zA-Z0-9_]*)", "m")
    all_files_there = true
    for match in eachmatch(pattern, content)
        # check if file already exists
        save_name = match.captures[1]
        if !isfile(joinpath(location, save_name))
            if haskey(ENV, "CI")
                throw(ErrorException("Running $fname in CI pipeline, missing $save_name at $location"))
            end
            all_files_there = false
            break
        end
    end
    # if all pre-rendered files exist, do nothing
    if all_files_there
        @info "Found pre-rendered images for $fname"
        return nothing
    end
    # if files missing, run script
    @warn "Detected missing pre-rendered file/s for script $fname, running..."
    include(fname)
    # scan for figs, copy to src folder
    for match in eachmatch(pattern, content)
        save_name = match.captures[1]
        path_name = dirname(cname)
        @info "Copying $save_name from \\build to \\src"
        cp(joinpath(path_name, save_name), joinpath(location, save_name); force=true)
    end
    return nothing
end

end