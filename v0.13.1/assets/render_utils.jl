function render_sphere!(ax, pos, r; kwargs...)
    # catch small radii render bug by capping min. radius
    if r < 1e-5
        r = 1e-5
    end
    render!(ax, BMO.SphereSDF(Point3{Float64}(pos), Float64(r)); transparency=true, kwargs...)
end

lens_color() = RGBf(0.678, 0.847, 0.902)
lens_color(alpha) = RGBAf(0.678, 0.847, 0.902, alpha)