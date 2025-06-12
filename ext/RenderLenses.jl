function render!(
        axis::_RenderEnv,
        css::BMO.ConcaveSphericalSurfaceSDF;
        # Makie kwargs
        color = :white,
        kwargs...
    )
    _radius = BMO.diameter(css)/2
    v = LinRange(0, 2π, 100)
    r = LinRange(1e-12, 1, 100) .^ (1/2) * _radius
    # Calculate beam surface at origin along y-axis, swap w and u
    y = sqrt.(BMO.radius(css)^2 .- r.^2) .- BMO.radius(css)
    u = y
    w = collect(r)
    # Close conture
    push!(u, 0) # push!(u, 0, 0)
    push!(w, _radius) # push!(w, radius, 1e-12)
    X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
    Y = [u for u in u, v in v]
    Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
    # Transform into global coords
    R = BMO.orientation(css)
    P = position(css)
    Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ P[1]
    Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ P[2]
    Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ P[3]
    surface!(axis, Xt, Yt, Zt; colormap = [color, color], kwargs...)
    return nothing
end

function render!(
        axis::_RenderEnv,
        css::BMO.ConvexSphericalSurfaceSDF;
        # Makie kwargs
        color = :white,
        kwargs...
    )
    v = LinRange(0, 2π, 100)
    r = LinRange(1e-12, 1, 100) .^ (1/2) * BMO.diameter(css)/2
    # Calculate beam surface at origin along y-axis, swap w and u
    y = -(sqrt.(BMO.radius(css)^2 .- r.^2) .- BMO.radius(css))
    u = y
    w = collect(r)
    # Close conture
    # push!(u, 0)
    # push!(w, 1e-12)
    X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
    Y = [u for u in u, v in v]
    Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
    # Transform into global coords
    R = BMO.orientation(css)
    P = position(css)
    Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ P[1]
    Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ P[2]
    Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ P[3]
    surface!(axis, Xt, Yt, Zt; colormap = [color, color], kwargs...)
    return nothing
end

function render!(
        axis::_RenderEnv,
        asp::BMO.AbstractAsphericalSurfaceSDF;
        color = :white,
        kwargs...
    )
    radius = asp.diameter / 2
    v = LinRange(0, 2π, 100)
    r = LinRange(1e-12, radius, 50)
    # Calculate beam surface at origin along y-axis, swap w and u
    y = BMO.aspheric_equation.(r, Ref(asp))
    u = y
    w = collect(r)
    if isa(asp, BMO.ConvexAsphericalSurfaceSDF)
        push!(u, u[end])
        push!(w, 1e-12)
    elseif isa(asp, BMO.ConcaveAsphericalSurfaceSDF)
        push!(u, 0, 0)
        push!(w, radius, 1e-12)
    else
        @warn "No suitable render fct. for $(typeof(asp))"
        return nothing
    end
    X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
    Y = [u for u in u, v in v]
    Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
    # Transform into global coords
    R = BMO.orientation(asp)
    P = position(asp)
    Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ P[1]
    Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ P[2]
    Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ P[3]
    surface!(axis, Xt, Yt, Zt; colormap = [color, color], kwargs...)
    return nothing
end