using BeamletOptics, GLMakie

function get_view(ls::LScene)
    cam = ls.scene.camera_controls
    eye = cam.eyeposition[]
    lookat = cam.lookat[]
    up = cam.upvector[]
    return eye, lookat, up
end

function set_view(ls::LScene, eye, lookat, up)
    cam = ls.scene.camera_controls
    cam.eyeposition[] = Vec3f(eye...)
    cam.lookat[] = Vec3f(lookat...)
    cam.upvector[] = Vec3f(up...)
    update_cam!(ls.scene, cam)
end

# Parameters
λ = 633e-9          # 633 nm (HeNe Laser)
a = 10e-6           # Slit width (10 µm)
d = 100e-6          # Slit separation (100 µm)
L = 0.2             # Propagation distance (0.2 m)
W = 0.5e-3          # Tighter window (0.5 mm)
nx = 251            # dx = 2 µm (10 beamlets per slit)

# 2. Define the Aperture Mask (Double Slit)
x = LinRange(-W / 2, W / 2, nx)
z = LinRange(-W / 2, W / 2, nx)
amplitude = zeros(nx, nx)

for i in 1:nx, j in 1:nx
    # Slits are long along z, narrow along x
    in_slit1 = abs(x[i] - d / 2) < a / 2
    in_slit2 = abs(x[i] + d / 2) < a / 2

    if (in_slit1 || in_slit2)
        amplitude[i, j] = 1.0
    end
end

phase = zeros(nx, nx)
dir = [0.0, 1.0, 0.0] # Propagation along Y
##
# Decompose into Astigmatic Beamlets
@info "Decomposing Double Slit into AGBs..."
# We explicitly define the basis to ensure Grid X = World X
e1 = [1.0, 0.0, 0.0]
e2 = [0.0, 0.0, 1.0]
beams = WavefrontBeamletDecomposition(x, z, amplitude, phase, dir, λ;
    threshold = 1e-3, overlap = 3.0, basis = (e1, e2))

@info "Number of beamlets generated: $(length(beams))"

# Define the System and Propagate
target_pd = Detector(30e-3)
translate_to3d!(target_pd, [0.0, L, 0.0])
system = System([target_pd])

@info "Propagating to screen at L = $L m..."
solve_system!(system, beams)

# Extract Results
@info "Calculating Intensity Pattern..."
n_eval = 500
# Evaluate over a 30mm window to see many fringes
xs_eval, zs_eval, I = intensity(target_pd; n = n_eval,
    x_min = -15e-3, x_max = 15e-3, z_min = -15e-3, z_max = 15e-3)

##
# Visualization
fig = Figure(size = (800, 600), fontsize = 18)
display(fig)

# Axis 1: 3D Scene (The Setup)
ax3d = LScene(fig[1:2, 1], show_axis = false)

# We create a mesh grid for the surface to ensure correct orientation in 3D
X_mask = [xv for xv in x, zv in z]
Z_mask = [zv for xv in x, zv in z]
Y_mask = fill(0.0, size(X_mask))
surface!(ax3d, X_mask, Y_mask, Z_mask,
    color = amplitude, colormap = :binary, transparency = true)

# Render Intensity Results (at y = L)
X_res = [xv for xv in xs_eval, zv in zs_eval]
Z_res = [zv for xv in xs_eval, zv in zs_eval]
Y_res = fill(L, size(X_res))
I_norm = I ./ (maximum(I) + 1e-12) # Avoid division by zero
plt_res = surface!(ax3d, X_res, Y_res, Z_res,
    color = I_norm, colormap = :magma,
    transparency = false,
    shading = NoShading,
    depth_shift = -1e-3) # Subtle shift to win Z-fighting without clipping

# Render a few representative beamlets
step = max(1, length(beams) ÷ 15)
for i in 1:step:length(beams)
    render!(ax3d, beams[i], color = (:cyan, 0.05), show_beams = false, flen = 0.0)
end

# Render the system components (Detector frame)
#render!(ax3d, system)

# Axis 2: 2D Intensity Map (Zoomed)
ax2d = Axis(fig[1, 2],
    title = "Diffraction Pattern (Detector Plane)",
    xlabel = "X [mm]", ylabel = "Z [mm]", yticks = -1:2:1,
    aspect = DataAspect())
hm = heatmap!(ax2d, xs_eval .* 1000, zs_eval .* 1000, I, colormap = :magma)
ylims!(-1.1, 1.1)

# Axis 3: 1D Cross-section
ax1d = Axis(fig[2, 2],
    title = "Interference Fringes (X-slice)",
    xlabel = "X Position [mm]", ylabel = "Intensity",
    yticksvisible = false, yticklabelsvisible = false)
lines!(ax1d, xs_eval .* 1000, I[:, size(I, 2) ÷ 2], color = :red, linewidth = 2)

# This makes the 1.2m distance look shorter so we can zoom in on X/Z details
scale!(ax3d.scene, 1.0, 0.15, 1.0)

# Set a nice default camera view
cl = ([0.005679703099338854, -0.010364784764641269, 0.002815435213044604],
    [-0.003708494535423113, 0.026029085430702627, 0.0009259786679933707],
    [-0.021598444746385476, 0.046277584411491254, 0.998695114799569])

rowsize!(fig.layout, 1, Relative(1 // 4))

set_view(ax3d, cl...)
save("agb_doubleslit_experiment.png", fig; px_per_unit = 8, update = false)
