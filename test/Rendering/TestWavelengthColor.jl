module TestWavelengthColor

using BeamletOptics
using Makie
using Test
using GeometryBasics: coordinates

const BMO = BeamletOptics

const Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

"""The color attribute of the plot `p` as a vector of colors."""
_colors(p) = to_value(p.color)

"""The `linesegments` plot of the handle `h`."""
_lines(h) = only(filter(p -> p isa Makie.LineSegments, BMO.render_plots(h)))

@testset "wavelength_color" begin
    r, g, b = wavelength_color(450e-9)
    @test b > r && b > g
    r, g, b = wavelength_color(550e-9)
    @test g > r && g > b
    r, g, b = wavelength_color(650e-9)
    @test r > g && r > b
    @test wavelength_color(650e-9) == (1.0, 0.0, 0.0)
    for λ in range(200e-9, 2000e-9; length = 200)
        @test all(c -> 0 <= c <= 1, wavelength_color(λ))
    end
    # dim, but not black, outside of the visible range, and continuous at its ends
    @test wavelength_color(300e-9) == wavelength_color(380e-9)
    @test wavelength_color(1064e-9) == wavelength_color(780e-9)
    @test 0 < maximum(wavelength_color(1064e-9)) < 0.5
    @test 0 < maximum(wavelength_color(300e-9)) < 0.5
    @test wavelength_color(400e-9) != wavelength_color(500e-9)
    @test wavelength_color(Float32(550e-9)) isa NTuple{3, Float64}
    @test_throws ArgumentError wavelength_color(0.0)
    @test_throws ArgumentError wavelength_color(-1e-9)
    @test_throws ArgumentError wavelength_color(NaN)
end

@testset "color = :wavelength" begin
    λs = [450e-9, 550e-9, 650e-9]
    rgba(λ, α = 1) = Makie.RGBAf(wavelength_color(λ)..., α)

    # a dispersive prism that splits the beams, a group of beams of different wavelengths
    prism = RightAnglePrism(30e-3, 20e-3, λ -> 1.5 + 4e-15 / λ^2)
    translate3d!(prism, [0, 20e-3, 0])
    # a beamsplitter whose child beams must keep the color of their wavelength
    bs = RoundThinBeamsplitter(25e-3; reflectance = 0.5)
    translate3d!(bs, [0, 2e-3, 0])
    zrotate3d!(bs, deg2rad(30))
    sys = System([bs, prism])
    beams = [Beam([0.0, 0, 0], [0.0, 1.0, 0.0], λ) for λ in λs]
    group = CollimatedSource(beams, 1e-3, [0.0, 0, 0], [0.0, 1.0, 0.0])
    solve_system!(sys, group)

    @testset "rays, beams and groups" begin
        ax = LScene(Figure()[1, 1])
        ray = Ray([0.0, 0, 0], [0.0, 1.0, 0.0], 450e-9)
        h = live_render!(ax, ray; color = :wavelength)
        @test _colors(_lines(h)) == [rgba(450e-9), rgba(450e-9)]

        h = live_render!(ax, beams[3]; color = :wavelength)
        pts = to_value(_lines(h).arg1)
        cols = _colors(_lines(h))
        @test length(cols) == length(pts)
        @test all(==(rgba(650e-9)), cols)

        # every beam of the group has its color, also after the prism
        h = live_render!(ax, group; color = :wavelength, render_every = 1)
        cols = _colors(_lines(h))
        @test length(cols) == length(to_value(_lines(h).arg1))
        @test Set(cols) == Set(rgba.(λs))
        @test all(b -> !isempty(b.children), beams)
        @test cols[1:2:end] == cols[2:2:end]
        nrays(b) = sum(length(BMO.rays(c)) for c in BMO.PreOrderDFS(b))
        @test count(==(rgba(550e-9)), cols) == 2 * nrays(beams[2])

        # opacity, show_pos, update after a change of the beam
        h = live_render!(ax, group; color = (:wavelength, 0.3), render_every = 1, show_pos = true)
        @test all(c -> c.alpha ≈ 0.3f0, _colors(_lines(h)))
        scat = only(filter(p -> p isa Makie.Scatter, BMO.render_plots(h)))
        @test length(_colors(scat)) == length(to_value(scat.arg1))
        update_render!(h)
        @test length(_colors(_lines(h))) == length(to_value(_lines(h).arg1))

        # render! draws the same
        n0 = length(ax.scene.plots)
        render!(ax, group; color = :wavelength, render_every = 1)
        @test Set(_colors(ax.scene.plots[n0 + 1])) == Set(rgba.(λs))
        # other colors are unchanged
        h = live_render!(ax, beams[1]; color = :red)
        @test _colors(_lines(h)) == Makie.to_color(:red)
    end

    @testset "Gaussian beamlets" begin
        ax = LScene(Figure()[1, 1])
        for λ in (450e-9, 650e-9)
            gauss = GaussianBeamlet([0.0, 0, 0], [0.0, 1.0, 0.0], λ, 1e-3)
            solve_system!(sys, gauss)
            h = live_render!(ax, gauss; color = (:wavelength, 0.5), r_res = 8, z_res = 6)
            mesh = only(filter(p -> p isa Makie.Mesh, BMO.render_plots(h)))
            cols = _colors(mesh)
            @test length(cols) == length(coordinates(to_value(mesh.arg1)))
            @test all(==(rgba(λ, 0.5)), cols)
            update_render!(h)
            @test length(_colors(mesh)) == length(coordinates(to_value(mesh.arg1)))
            render!(ax, gauss; color = :wavelength, r_res = 8, z_res = 6, show_beams = true)
        end

        agb = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 550e-9, 1e-3; support = [0, 0, 1])
        solve_system!(sys, agb; check_invariant = false)
        h = live_render!(ax, agb; color = :wavelength, r_res = 8, z_res = 6, show_waist = true)
        mesh = only(filter(p -> p isa Makie.Mesh, BMO.render_plots(h)))
        @test length(_colors(mesh)) == length(coordinates(to_value(mesh.arg1)))
        @test all(==(rgba(550e-9)), _colors(mesh))

        bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 4e-3, 450e-9, 1e-3; n_grid = 3)
        solve_system!(sys, bg; check_invariant = false)
        h = live_render!(ax, bg; color = :wavelength, render_every = 2)
        mesh = only(filter(p -> p isa Makie.Mesh, BMO.render_plots(h)))
        @test length(_colors(mesh)) == length(coordinates(to_value(mesh.arg1)))
        @test all(==(rgba(450e-9)), _colors(mesh))
    end
end

end # module
