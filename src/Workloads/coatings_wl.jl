@setup_workload begin
    λ = 1000e-9
    @compile_workload begin
        nc = sqrt(1.5)
        d = λ / (4 * nc)
        coat_thin = ThinFilmCoating(nc, d; behavior = Transmissive())
        
        n_lens = λ -> 1.5
        base_lens = SphericalLens(100e-3, -100e-3, 5e-3, 10e-3, n_lens)
        coated_lens = base_lens |> with_coatings(front = coat_thin)

        system = System([coated_lens])

        # Trace plain Ray
        ray_lr = Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ)
        beam_lr = Beam(ray_lr)
        solve_system!(system, beam_lr; retrace = false)

        # Trace PolarizedRay
        pray = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])
        pbeam = Beam(pray)
        solve_system!(system, pbeam; retrace = false)

        # Trace GaussianBeamlet
        gauss = GaussianBeamlet([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, 1e-3)
        solve_system!(system, gauss; retrace = false)
    end
end
