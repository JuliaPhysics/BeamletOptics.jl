# Run workloads on user machine only
if get(ENV, "CI", "false") == "false"
    let 
        include("michelson_wl.jl")
        include("asphere_wl.jl")
    end
end