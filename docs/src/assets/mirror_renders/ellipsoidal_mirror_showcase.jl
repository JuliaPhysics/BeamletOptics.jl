using GLMakie, BeamletOptics

GLMakie.activate!(; ssao = true)

const BMO = BeamletOptics
const mm = 1e-3

## Ellipsoidal mirror imaging its far focus into its near focus
# semi-major axis a and focal distance c of the parent ellipse, foci at y = -(a - c) and y = -(a + c)
a_val = 100mm
c_val = 60mm
D_val = 100mm

s_near = a_val - c_val
s_far = a_val + c_val

# vertex at the origin, opens towards -y
em = EllipsoidalMirror(s_near, s_far, D_val; thickness = 5mm)

F_near = [0, -s_near, 0]
F_far = [0, -s_far, 0]

system = StaticSystem([em])

## Fan of rays from the far focus, aimed at points across the mirror aperture (x-y plane)
R_val = BMO.radius(BMO.shape(em))
k_val = BMO.shape(em).k

beams = Beam{Float64, Ray{Float64}}[]
for z in range(-0.8, 0.8; length = 20) .* D_val / 2
    target = [0, -BMO._conic_sag(abs(z), R_val, k_val), z]
    local beam = Beam(Ray(F_far, target - F_far))
    solve_system!(system, beam)
    push!(beams, beam)
end

##
cview = [
 0.0896927   0.995969    -1.0365e-16   0.0613752
 0.0331685  -0.00298701   0.999445     0.000141382
 0.995417   -0.0896429   -0.0333027   -0.149613
 0.0         0.0          0.0          1.0
]

fig = Figure(size = (600, 350))
display(fig)
ax = LScene(fig[1, 1], show_axis = false)

render!(ax, em; transparency = true, alpha = 0.2)
for beam in beams
    # the reflected rays are untruncated, let them run through the near focus a bit
    render!(ax, beam; color=:blue, alpha=0.5, show_pos = true, flen = 1.5s_near)
end
scatter!(ax, [Point3f(F_near), Point3f(F_far)]; color = :red, markersize = 12)

set_view(ax, cview)
save("ellipsoidal_mirror_showcase.png", fig; px_per_unit = 4, update = false)
