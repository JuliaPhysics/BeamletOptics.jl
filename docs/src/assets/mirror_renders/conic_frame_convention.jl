using GLMakie, BeamletOptics

GLMakie.activate!()

const BMO = BeamletOptics

## Frame-convention sketch of a ConicSDF segment (2D cross-section in the x-y plane)
R_sketch = 0.2
k_sketch = -0.5
x_off_sketch = 0.06
D_sketch = 0.08

# surface in the segment frame: y = -(Z(|x + x_off|) - Z(x_off)), segment centre at the origin
Z_off = BMO._conic_sag(x_off_sketch, R_sketch, k_sketch)
surf_y(x) = -(BMO._conic_sag(abs(x + x_off_sketch), R_sketch, k_sketch) - Z_off)

xs = range(-D_sketch / 2, D_sketch / 2; length = 200)
parent_xs = range(-x_off_sketch - 0.09, D_sketch / 2 + 0.02; length = 300)

# parent vertex and near focus on the parent axis x = -x_off, focal distance R/(1 + e) with e = sqrt(-k)
vertex_pt = Point2f(-x_off_sketch, Z_off)
focus_pt = Point2f(-x_off_sketch, Z_off - R_sketch / (1 + sqrt(-k_sketch)))

fig2 = Figure(size = (600, 450))
display(fig2)
ax3 = Axis(fig2[1, 1]; xlabel = "x [m]", ylabel = "y [m]", aspect = DataAspect())

# parent axis
vlines!(ax3, -x_off_sketch; color = :grey, linestyle = :dashdot, label = "parent axis")
lines!(ax3, parent_xs, surf_y.(parent_xs); color = :grey, linestyle = :dash, label = "parent conic")
lines!(ax3, xs, surf_y.(xs); color = :black, linewidth = 3, label = "segment aperture")

# origin = segment centre on the surface, with the -y opening direction
scatter!(ax3, Point2f(0, 0); color = :blue, markersize = 12)
text!(ax3, 0.0, 0.0; text = " origin\n (segment centre)", align = (:left, :bottom))
arrows2d!(ax3, [Point2f(0, 0)], [Vec2f(0, -0.04)]; color = :blue)
text!(ax3, 0.0, -0.03; text = " opens towards -y", color = :blue, align = (:left, :center))

# parent vertex and near focus
scatter!(ax3, vertex_pt; color = :red, markersize = 12)
text!(ax3, vertex_pt; text = "parent vertex\n(-x_off, +Z(x_off)) ", align = (:right, :bottom))
scatter!(ax3, focus_pt; color = :green, markersize = 12)
text!(ax3, focus_pt; text = " near focus", align = (:left, :center))

# off-axis distance between the parent axis and the segment centre
y_dim = 0.035
bracket!(ax3, -x_off_sketch, y_dim, 0.0, y_dim; text = "x_off", orientation = :up, color = :black)

ylims!(ax3, focus_pt[2] - 0.02, 0.08)
axislegend(ax3; position = :rb)

save("conic_frame_convention.png", fig2; px_per_unit = 8, update = false)
