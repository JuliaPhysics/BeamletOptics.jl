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
        # a single color is passed on by render!, a handle draws it per vertex, see render_settings!
        n0 = length(ax.scene.plots)
        render!(ax, beams[1]; color = :red)
        @test _colors(ax.scene.plots[n0 + 1]) == Makie.to_color(:red)
        h = live_render!(ax, beams[1]; color = :red)
        @test length(_colors(_lines(h))) == length(to_value(_lines(h).arg1))
        @test all(==(Makie.RGBAf(Makie.to_color(:red))), _colors(_lines(h)))
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

"""The number of vertices of the positions (`arg1`) of the plot `p`: of a mesh, points or segments."""
_nvertices(p) = (a = to_value(p.arg1); a isa AbstractVector ? length(a) : length(coordinates(a)))

"""
    _transient_lengths(h, solve)

Calls `solve()` and `update_render!(h)` and returns the pairs `(vertices, colors)` of the `color = :wavelength`
plots of `h` that Makie sees whenever the positions or the colors of one of them are notified, i.e. at
every step of the update.
"""
function _transient_lengths(h, solve)
    seen = Tuple{Int, Int}[]
    plots = filter(p -> to_value(p.color) isa AbstractVector{<:Makie.RGBAf}, BMO.render_plots(h))
    for p in plots
        check(_) = push!(seen, (_nvertices(p), length(to_value(p.color))))
        on(check, p.arg1)
        on(check, p.color)
    end
    solve()
    update_render!(h)
    return seen
end

@testset "color = :wavelength: positions and colors are updated together" begin
    # a mirror that is moved in and out of the beam changes the number of segments on every solve
    lens = SphericalLens(50e-3, -50e-3, 10e-3, 25e-3, 1.5168)
    ax = LScene(Figure()[1, 1])

    @testset "$name" for (name, make, kwargs) in (
            ("Beam", () -> Beam([0.0, -20e-3, 0.0], [0.0, 1.0, 0.0], 550e-9), (; show_pos = true, render_every = 1)),
            ("GaussianBeamlet", () -> GaussianBeamlet([0.0, -20e-3, 0.0], [0.0, 1.0, 0.0], 550e-9, 1e-3),
                (; r_res = 8, z_res = 6)),
            ("AstigmaticGaussianBeamlet",
                () -> AstigmaticGaussianBeamlet([0, -20e-3, 0], [0, 1, 0], 550e-9, 1e-3; support = [0, 0, 1]),
                (; r_res = 8, z_res = 6, show_waist = true)),
        )
        mir = RoundPlanoMirror(25.4e-3, 5e-3)
        translate3d!(mir, [0, 60e-3, 0])
        zrotate3d!(mir, deg2rad(45))
        sys = System([lens, mir])
        thing = make()
        solve_system!(sys, thing; check_invariant = false)
        h = live_render!(ax, thing; color = :wavelength, kwargs...)
        sizes = Int[]
        for d in (10.0, -10.0, 10.0, -10.0)
            seen = _transient_lengths(h, () -> begin
                translate3d!(mir, [d, 0, 0])
                solve_system!(sys, thing; check_invariant = false)
            end)
            push!(sizes, maximum(first, seen; init = 0))
            @test !isempty(seen)
            @test all(l -> l[1] == l[2], seen)
        end
        # the vertex count did change
        @test length(unique(sizes)) > 1
    end
end

@testset "render_settings!: color" begin
    rgba(λ, α = 1) = Makie.RGBAf(wavelength_color(λ)..., α)
    red = Makie.RGBAf(Makie.to_color(:red))
    lens = SphericalLens(50e-3, -50e-3, 10e-3, 25e-3, 1.5168)
    sys = System([lens])
    ax = LScene(Figure()[1, 1])
    # positions and colors of every colored plot of `h` have the same length
    consistent(h) = all(p -> length(_colors(p)) == _nvertices(p),
        filter(p -> _colors(p) isa AbstractVector, BMO.render_plots(h)))

    @testset "rays, beams and groups" begin
        beams = [Beam([x, -20e-3, 0.0], [0.0, 1.0, 0.0], λ) for (x, λ) in ((-1e-3, 450e-9), (1e-3, 650e-9))]
        group = CollimatedSource(beams, 1e-3, [0.0, -20e-3, 0], [0.0, 1.0, 0.0])
        solve_system!(sys, group)
        h = live_render!(ax, group; color = :wavelength, render_every = 1, show_pos = true)
        plots = copy(BMO.render_plots(h))
        nplots = length(ax.scene.plots)
        @test Set(_colors(_lines(h))) == Set(rgba.((450e-9, 650e-9)))

        # a single color, without new plots
        @test BMO.render_settings!(h; color = :red) === h
        @test BMO.render_settings(h).color === :red
        @test all(p -> all(==(red), _colors(p)), BMO.render_plots(h))
        @test consistent(h)
        @test length(BMO.render_plots(h)) == length(plots)
        @test all(a === b for (a, b) in zip(BMO.render_plots(h), plots))
        @test length(ax.scene.plots) == nplots
        # it is kept by the updates
        update_render!(h)
        BMO.render_settings!(h; flen = 0.2)
        @test all(p -> all(==(red), _colors(p)), BMO.render_plots(h))
        @test consistent(h)

        # with an opacity, and back to the wavelengths
        BMO.render_settings!(h; color = (:green, 0.25))
        @test all(c -> c.alpha ≈ 0.25f0, _colors(_lines(h)))
        BMO.render_settings!(h; color = (:wavelength, 0.5))
        @test BMO.render_settings(h).color == (:wavelength, 0.5)
        @test Set(_colors(_lines(h))) == Set(rgba.((450e-9, 650e-9), 0.5))
        @test all(a === b for (a, b) in zip(BMO.render_plots(h), plots))
        @test length(ax.scene.plots) == nplots

        # the colors follow the other settings
        BMO.render_settings!(h; render_every = 2)
        @test Set(_colors(_lines(h))) == Set((rgba(450e-9, 0.5),))
        @test consistent(h)

        # a handle of a single color takes the wavelengths
        h = live_render!(ax, group; render_every = 1)
        @test BMO.render_settings(h).color === :blue
        BMO.render_settings!(h; color = :wavelength)
        @test Set(_colors(_lines(h))) == Set(rgba.((450e-9, 650e-9)))
        @test consistent(h)
    end

    @testset "Gaussian beamlets" begin
        gauss = GaussianBeamlet([0.0, -20e-3, 0.0], [0.0, 1.0, 0.0], 650e-9, 1e-3)
        solve_system!(sys, gauss)
        h = live_render!(ax, gauss; r_res = 8, z_res = 6, show_beams = true)
        mesh = only(filter(p -> p isa Makie.Mesh, BMO.render_plots(h)))
        gen = filter(p -> p isa Makie.LineSegments, BMO.render_plots(h))
        gencolors = [to_value(p.color) for p in gen]
        @test all(==(red), _colors(mesh))
        BMO.render_settings!(h; color = (:wavelength, 0.5), z_res = 4)
        @test all(==(rgba(650e-9, 0.5)), _colors(mesh))
        @test consistent(h)
        BMO.render_settings!(h; color = :green)
        @test all(==(Makie.RGBAf(Makie.to_color(:green))), _colors(mesh))
        # the overlay of show_beams keeps its colors
        @test [to_value(p.color) for p in gen] == gencolors

        bg = CollimatedGaussianBeamletSource([0, -20e-3, 0], [0, 1, 0], 4e-3, 450e-9, 1e-3; n_grid = 3)
        solve_system!(sys, bg; check_invariant = false)
        h = live_render!(ax, bg; render_every = 2, show_waist = true)
        BMO.render_settings!(h; color = :wavelength, render_every = 1)
        @test all(p -> all(==(rgba(450e-9)), _colors(p)),
            filter(p -> p isa Union{Makie.Mesh, Makie.Scatter}, BMO.render_plots(h)))
        @test consistent(h)
    end

    @testset "errors" begin
        beam = Beam([0.0, -20e-3, 0.0], [0.0, 1.0, 0.0], 550e-9)
        solve_system!(sys, beam)
        h = live_render!(ax, beam)
        for color in (1, :nocolor, [:red, :blue], (:red, :blue))
            @test_throws ArgumentError BMO.render_settings!(h; color)
        end
        @test_throws ArgumentError BMO.render_settings!(h; color = :red, flen = -1.0)
        @test BMO.render_settings(h).color === :blue
        # colors that are passed on to Makie are not managed by the handle
        n = length(to_value(_lines(h).arg1))
        hv = live_render!(ax, beam; color = fill(Makie.RGBAf(0, 1, 0, 1), n))
        @test_throws ArgumentError BMO.render_settings!(hv; color = :red)
        BMO.render_settings!(hv; flen = 0.2)
    end
end

end # module
