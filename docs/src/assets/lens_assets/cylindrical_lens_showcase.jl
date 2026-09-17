using GLMakie, BeamletOptics

##
r = 5.2e-3  # radius
d = 10e-3   # diameter/width of the cylindric portion
h = 20e-3   # height/length of the cylinder
ct = 5.9e-3 # center thickness
cylindrical_lens = Lens(
    CylindricalSurface(r, d, h),
    ct,
    n -> 1.517
)

fig_cylindrical = Figure(size=(600, 450))
ax = Axis3(fig_cylindrical[1,1], aspect=:data, azimuth=-pi/4, elevation=deg2rad(30))
hidedecorations!(ax)
hidespines!(ax)
render!(ax, cylindrical_lens)

save("cylindrical_lens_showcase.png", fig_cylindrical)

##
radius = -15.538e-3
diameter = 25e-3
height = 50e-3
conic_constant = -1.0

acylindrical_lens = Lens(
    BeamletOptics.AcylindricalSurface(
        radius,
        diameter,
        height,
        conic_constant,
        [0, 1.1926075e-5*(1e3)^3, -2.9323497e-9*(1e3)^5, -1.8718889e-11*(1e3)^7, -1.7009961e-14*(1e3)^9, 3.5481542e-17*(1e3)^11, 6.5241296e-20*(1e3)^13]
    ),
    7.5e-3,
    n -> 1.777
)

fig_acylindrical = Figure(size=(600, 450))
ax = Axis3(fig_acylindrical[1,1], aspect=:data, azimuth=-pi/4, elevation=deg2rad(30))
hidedecorations!(ax)
hidespines!(ax)
render!(ax, acylindrical_lens)

save("acylindrical_lens_showcase.png", fig_acylindrical)
