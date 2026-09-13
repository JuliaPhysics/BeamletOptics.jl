"""
    ReflectiveGrating{T} <: BeamletOptics.AbstractObject{T}

A planar reflective diffraction grating, modelled as a [`BeamletOptics.MultiShape`](@ref)
object with two parts:

- `grating`: a flat `Mesh` representing the ruled front face, where diffraction occurs
- `substrate`: an opaque `BoxSDF` behind the front face, which terminates any ray that
  reaches it

Splitting the component this way is what makes the grating one-sided: rays arriving from
the back are absorbed by the substrate instead of being diffracted through it.

# Fields

- `grating`: front face mesh
- `substrate`: opaque backing volume
- `groove_density`: grooves per meter, i.e. `1/d`
- `order`: diffraction order `m`
"""
struct ReflectiveGrating{T} <: BMO.AbstractObject{T}
    grating::BMO.Mesh{T}
    substrate::BMO.BoxSDF{T}
    groove_density::T           # lines/meter
    order::Int8
end

BMO.shape_trait_of(::ReflectiveGrating) = BMO.MultiShape()

BMO.shape(rg::ReflectiveGrating) = (rg.grating, rg.substrate)

"""
    RectangularReflectiveGrating(width, height, thickness, groove_density, order)

Spawns a rectangular [`ReflectiveGrating`](@ref) at the global origin. The ruled face lies
in the local xz-plane, the substrate extends along the local +y-axis and the local x-axis
is the dispersion direction.

# Arguments

- `width`: extent along the local x-axis (dispersion direction) in [m]
- `height`: extent along the local z-axis in [m]
- `thickness`: substrate thickness in [m]
- `groove_density`: grooves per meter, e.g. `1.2e6` for a 1200 lines/mm grating
- `order`: diffraction order `m`
"""
function RectangularReflectiveGrating(
        width::Real,
        height::Real,
        thickness::Real,
        groove_density::Real,
        order::Int
    )
    grating = BMO.RectangularFlatMesh(width, height)
    substrate = BMO.BoxSDF(width, thickness, height)
    # offset the substrate by a hair so that the two shapes do not share a surface
    translate3d!(substrate, [0, thickness/2 + 1e-6, 0])
    return ReflectiveGrating(grating, substrate, groove_density, Int8(order))
end

function BMO.render!(ax::LScene, gr::ReflectiveGrating; kwargs...)
    BMO.render!(ax, gr.grating; color=:orange, kwargs...)
    BMO.render!(ax, gr.substrate; color=:white, transparency=true, kwargs...)
    return nothing
end

function BMO.interact3d(
        ::BMO.AbstractSystem,
        gr::ReflectiveGrating,
        ::Beam{T, R},
        ray::R
    ) where {T <: Real, R <: Ray{T}}

    # if substrate hit, stop trace
    if BMO.shape(BMO.intersection(ray)) === gr.substrate
        return nothing
    end

    normal = BMO.normal3d(BMO.intersection(ray))
    v_in = BMO.direction(ray)
    lambda = BMO.wavelength(ray)

    dispersion_dir = BMO.orientation(gr)[:,1]

    # Project dispersion vector onto the local tangent plane of the grating surface
    g_proj = dispersion_dir - dot(dispersion_dir, normal) * normal
    g_hat = normalize(g_proj)

    # Decompose incident direction into tangential and normal components
    v_in_dot_n = dot(v_in, normal)
    v_in_normal = v_in_dot_n * normal
    v_in_tangent = v_in - v_in_normal

    # Vector Grating Equation: v_out_tangent = v_in_tangent + m * (λ / d) * g_hat
    v_out_tangent = v_in_tangent + (gr.order * lambda * gr.groove_density) * g_hat

    # Check for evanescent wave (ray diffracts past 90 degrees)
    tangent_sq = dot(v_out_tangent, v_out_tangent)
    if tangent_sq > one(T)
        @debug "ray misses grating mode / evanescent cutoff"
        return nothing
    end

    # Calculate reflected normal component (rebound opposite to incident sign)
    v_out_normal_mag = sqrt(max(zero(T), one(T) - tangent_sq))
    v_out_normal = (v_in_dot_n < 0 ? normal : -normal) * v_out_normal_mag

    # Reconstruct final ray direction and position
    ndir = normalize(v_out_tangent + v_out_normal)
    npos = BMO.position(ray) + BMO.length(ray) * BMO.direction(ray)

    return BMO.BeamInteraction{T, R}(
        nothing, # No secondary surface interaction hint required for single reflection
        BMO.Ray{T}(npos, ndir, nothing, lambda, BMO.refractive_index(ray))
    )
end
