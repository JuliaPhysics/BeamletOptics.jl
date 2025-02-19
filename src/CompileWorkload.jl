@setup_workload begin
    # setup dummy workload
    cm = 1e-2
    splitter_origin = [18.81cm, 23.5cm, 0]

    @compile_workload begin
        # execute workload
        NBK7 = SCDI.DiscreteRefractiveIndex([632.8e-9], [1.51509])

        # Mirror
        rpm = SCDI.RightAnglePrismMirror(25e-3, 25e-3)
        SCDI.zrotate3d!(rpm, deg2rad(45))
        SCDI.translate3d!(rpm, [0, 33.5cm, 0])

        mirror_assembly = SCDI.ObjectGroup([rpm])

        # Beamsplitter
        cbs = SCDI.CubeBeamsplitter(SCDI.inch, NBK7)
        SCDI.zrotate3d!(cbs, deg2rad(-90))

        splitter_assembly = SCDI.ObjectGroup([cbs])

        # Arms
        m1 = SCDI.RoundPlanoMirror(SCDI.inch, 5e-3)
        SCDI.zrotate3d!(m1, deg2rad(-90))
        SCDI.translate3d!(m1, [22cm, 0, 0])
        m2 = SCDI.RoundPlanoMirror(SCDI.inch, 5e-3)
        SCDI.zrotate3d!(m2, deg2rad(-90))
        SCDI.translate3d!(m2, [12cm, 0, 0])

        arm_1 = SCDI.ObjectGroup([m1])
        arm_2 = SCDI.ObjectGroup([m2])

        # PD
        pd = SCDI.Photodetector(8e-3, 200)
        SCDI.translate3d!(pd, [0, -12cm, 0])

        pd_assembly = SCDI.ObjectGroup([pd])

        system = SCDI.System(
            [
            mirror_assembly,
            splitter_assembly,
            arm_1,
            arm_2,
            pd_assembly
        ]
        )

        ##
        SCDI.translate_to3d!(mirror_assembly, [0, -10cm, 0])
        SCDI.translate_to3d!(splitter_assembly, [18.81cm, 23.5cm, 0])

        SCDI.translate_to3d!(arm_1, splitter_origin)
        SCDI.translate_to3d!(arm_2, splitter_origin)
        SCDI.translate3d!(arm_1, [3.81cm / 2, 0, 0])
        SCDI.translate3d!(arm_2, [0, 3.81cm / 2, 0])
        SCDI.zrotate3d!(arm_2, deg2rad(90))

        SCDI.translate_to3d!(pd_assembly, splitter_origin)
        SCDI.translate3d!(pd_assembly, [0, -3.81cm / 2, 0])

        ##
        beam = SCDI.GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 632.8e-9, 5e-4, M2 = 2)
        SCDI.solve_system!(system, beam)

        ##
        SCDI.reset_detector!(pd)
        SCDI.translate3d!(m1, [5e-9, 0, 0])
        SCDI.solve_system!(system, beam)
        p = SCDI.optical_power(pd)
    end

    @compile_workload begin
        m1 = SCDI.SquarePlanoMirror(SCDI.inch, 10e-3)
        m2 = SCDI.SquarePlanoMirror(SCDI.inch, 10e-3)
        bs1 = SCDI.ThinBeamsplitter(SCDI.inch, SCDI.inch)
        bs2 = SCDI.ThinBeamsplitter(SCDI.inch, SCDI.inch)
        pd = SCDI.Photodetector(SCDI.inch, 200)

        SCDI.translate3d!(bs1, [10cm, 0, 0])
        SCDI.translate3d!(m1, [20cm, 0, 0])
        SCDI.translate3d!(m2, [10cm, 10cm, 0])
        SCDI.translate3d!(bs2, [20cm, 10cm, 0])
        SCDI.translate3d!(pd, [30cm, 10cm, 0])

        SCDI.zrotate3d!(bs1, deg2rad(45))
        SCDI.zrotate3d!(bs2, deg2rad(45))
        SCDI.zrotate3d!(m1, deg2rad(45 + 180))
        SCDI.zrotate3d!(m2, deg2rad(45))
        SCDI.zrotate3d!(pd, deg2rad(90))

        system = SCDI.System([m1, m2, bs1, bs2, pd])

        beam = SCDI.GaussianBeamlet([0, 0, 0], [1, 0, 0], 632e-9, 2e-3)
        SCDI.reset_detector!(pd)
        SCDI.solve_system!(system, beam)
        p = SCDI.optical_power(pd)
    end

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
    n = 1.5

    @compile_workload begin
        NBK7 = SCDI.DiscreteRefractiveIndex([1000e-9], [1.5])
        NLAK22 = SCDI.DiscreteRefractiveIndex([1000e-9], [1.6])
        NSF10 = SCDI.DiscreteRefractiveIndex([1000e-9], [1.7])

        LB1811 = SCDI.SphericalLens(34.9e-3, -34.9e-3, 6.8e-3, SCDI.inch, NBK7)
        AC254_150_AB = SCDI.SphericalDoubletLens(
            87.9e-3, 105.6e-3, 1000, 6e-3, 3e-3, SCDI.inch, NLAK22, NSF10)

        aspherical_lens = SCDI.Lens(
            SCDI.generalized_lens_shape_constructor(R, Inf, ct, d;
                front_kind = :aspherical, front_k = k, front_coeffs = A
            ),
            _n -> n
        )

        SCDI.translate3d!(LB1811, [0, 10cm, 0])
        SCDI.translate3d!(AC254_150_AB, [0, 20cm, 0])
        SCDI.translate3d!(aspherical_lens, [0, 40cm, 0])

        system = SCDI.System([LB1811, AC254_150_AB, aspherical_lens])

        beam = SCDI.GaussianBeamlet([0, 0, 0], [0, 1, 0], 1000e-9, 2e-3)
    end
end