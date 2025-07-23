@setup_workload begin
    cm = 1e-2
    splitter_origin = [18.81cm, 23.5cm, 0]

    @compile_workload begin
        NBK7 = DiscreteRefractiveIndex([632.8e-9], [1.51509])

        # Mirror
        rpm = RightAnglePrismMirror(25e-3, 25e-3)
        zrotate3d!(rpm, deg2rad(45))
        translate3d!(rpm, [0, 33.5cm, 0])

        mirror_assembly = ObjectGroup([rpm])

        # Beamsplitter
        cbs = CubeBeamsplitter(inch, NBK7)
        zrotate3d!(cbs, deg2rad(-90))

        splitter_assembly = ObjectGroup([cbs])

        # Arms
        m1 = RoundPlanoMirror(inch, 5e-3)
        zrotate3d!(m1, deg2rad(-90))
        translate3d!(m1, [22cm, 0, 0])
        m2 = RoundPlanoMirror(inch, 5e-3)
        zrotate3d!(m2, deg2rad(-90))
        translate3d!(m2, [12cm, 0, 0])

        arm_1 = ObjectGroup([m1])
        arm_2 = ObjectGroup([m2])

        # PD
        pd = Photodetector(8e-3, 200)
        translate3d!(pd, [0, -12cm, 0])

        pd_assembly = ObjectGroup([pd])

        system = System([
            mirror_assembly,
            splitter_assembly,
            arm_1,
            arm_2,
            pd_assembly
        ])

        ##
        translate_to3d!(mirror_assembly, [0, -10cm, 0])
        translate_to3d!(splitter_assembly, [18.81cm, 23.5cm, 0])

        translate_to3d!(arm_1, splitter_origin)
        translate_to3d!(arm_2, splitter_origin)
        translate3d!(arm_1, [3.81cm / 2, 0, 0])
        translate3d!(arm_2, [0, 3.81cm / 2, 0])
        zrotate3d!(arm_2, deg2rad(90))

        translate_to3d!(pd_assembly, splitter_origin)
        translate3d!(pd_assembly, [0, -3.81cm / 2, 0])

        ##
        beam = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 632.8e-9, 5e-4, M2 = 2)
        solve_system!(system, beam)

        ##
        empty!(pd)
        translate3d!(m1, [5e-9, 0, 0])
        solve_system!(system, beam)
        p = optical_power(pd)
    end
end