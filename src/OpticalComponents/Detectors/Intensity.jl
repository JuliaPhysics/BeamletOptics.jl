electric_field(d::Detector; kwargs...) = electric_field(d, hits(d); kwargs...)

function electric_field(pd::Detector, hits::Vector{GaussianBeamletHit{G}}; n::Int=100, kwargs...) where G
    # Preallocate transforms
    T = transpose(orientation(shape(pd)))
    p = position(shape(pd))
    # Preallocate x-/y-buffers and field
    xs = ys = LinRange(-pd.edgln/2, pd.edgln/2, n)
    field = zeros(Complex{G}, n, n)
    # Calculate field superposition
    Threads.@threads for j in eachindex(ys) # FIXME row column major order?
        y = ys[j]
        @inbounds for i in eachindex(xs)
            x = xs[i] 
            # Transform point p on PD into world coordinates
            p1 = Point3(
                T[1, 1] * x + T[1, 3] * y + p[1],
                T[2, 1] * x + T[2, 3] * y + p[2],
                T[3, 1] * x + T[3, 3] * y + p[3]
            )
            # Add E-field contribution for each hit
            for hit in hits
                p0 = position(hit)
                d0 = direction(hit)
                l0 = length(hit)
                proj = projection_factor(hit)
                # Find projection of p1 onto Gaussian optical axis, i.e. local r and z
                l1 = dot(p1 - p0, d0)
                p2 = p0 + l1 * d0
                r = norm(p1 - p2)
                z = l0 + l1
                # Add field contribution, projection factor accounts for beam spot stretching
                field[i, j] += electric_field(hit.gauss, r, z) * sqrt(proj)
            end
        end
    end
    return xs, ys, field
end

function intensity(d::Detector; kwargs...)
    x, y, E = electric_field(d; kwargs...)
    return x, y, intensity.(E)
end