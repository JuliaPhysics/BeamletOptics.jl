@setup_workload begin
    ## Parameters for aspherical lens
    # radius
    R = 50.3583e-3
    # conic constant
    k = -0.789119
    # even aspheric coefficients up to 8th order
    A = [0, 2.10405e-7 * (1e3)^3, 1.76468e-11 * (1e3)^5, 1.02641e-15 * (1e3)^7]
    # center thickness
    ct = 10.2e-3
    # diameter
    d = 50e-3
    # refractive index of BK-7 @ 1310 nm (design wavelength)
    n = 1.5036

    @compile_workload begin
        lens = Lens(
            BeamletOptics.generalized_lens_shape_constructor(R, Inf, ct, d;
                front_kind = :aspherical, front_k = k, front_coeffs = A
            ),
            _n -> n
        )

        system = System(lens)

        for z in -0.02:0.001:0.02
            pos = [0.0, -0.05, z]
            dir = [0.0, 1.0, 0]
            ray = Ray(pos, dir)
            beam = Beam(ray)
            solve_system!(system, beam, r_max = 40)
        end
    end
end