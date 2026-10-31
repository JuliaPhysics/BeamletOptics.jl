module TestRenderMirrors

using BeamletOptics
using Makie
using GeometryBasics
using Test

const BMO = BeamletOptics

"""Returns the new plots of `render!(ax, x; kwargs...)`."""
function rendered_plots(x; kwargs...)
    ax = LScene(Figure()[1, 1])
    n0 = length(ax.scene.plots)
    render!(ax, x; kwargs...)
    return ax.scene.plots[(n0 + 1):end]
end

"""Pierced conic mirror whose bore is tilted about x by `θ` about its own position."""
function pierced_mirror(θ)
    mirror = ConicMirror(200e-3, -0.5, 50e-3; thickness = 6e-3, hole_diameter = 5e-3)
    tool = BMO.operands(BMO.shape(mirror))[2]
    p0 = copy(position(tool))
    xrotate3d!(tool, θ)
    @assert position(tool) == p0
    return mirror
end

"""Vertex coordinates of the mesh that `render!` draws for `mirror`."""
mesh_coordinates(mirror) = collect(GeometryBasics.coordinates(rendered_plots(mirror)[1][1][]))

@testset "Coaxiality of the analytic bore mesh" begin
    coaxial = pierced_mirror(0.0)
    @test BMO.shape(coaxial) isa BMO.DifferenceSDF{Float64, <:BMO.ConicSDF}
    reference = mesh_coordinates(coaxial)
    # within the tolerance of about 2.56°, parallel and anti-parallel bores get the analytic mesh
    @test isequal(mesh_coordinates(pierced_mirror(0.04)), reference)
    @test isequal(mesh_coordinates(pierced_mirror(π - 0.04)), reference)
    # beyond it, the bore is cut from the mesh of the conic
    @test !isequal(mesh_coordinates(pierced_mirror(0.05)), reference)
end

end
