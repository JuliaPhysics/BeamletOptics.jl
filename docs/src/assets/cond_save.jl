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
        pattern = Regex("save\\(\\s*\\\"([^\\\"]+)\\\"\\s*,\\s*([a-zA-Z_][a-zA-Z0-9_]*)")
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
            Box(fig[1, 1], color = :gray)
            save(save_name, fig)
            @info "Saving placeholder for $save_name (Size: $_width x $_height, fig_var_name=$fig_name)"
        end
    else
        @info "Running script for $fname"
        include(fname)
    end
end