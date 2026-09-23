using GLMakie, BeamletOptics

const BMO = BeamletOptics
const mm = 1e-3

##
s1 = BMO.SphereSDF(20mm)
s2 = BMO.CylinderSDF(5mm, 30mm)

# translate_to3d!(s1, [0, -5mm, 0])
# translate_to3d!(s2, [0, +5mm, 0])

##
cview = [
 -0.839316   0.543644  -1.94289e-16  -0.000655699
 -0.172994  -0.267081   0.948019      0.0014867
  0.515385   0.795688   0.318213     -0.0693439
  0.0        0.0        0.0           1.0
]

fig = Figure(size=(600, 600))
ax = LScene(fig[1,1], show_axis=false)
res = LScene(fig[1,2], show_axis=false)
tes = LScene(fig[2,1], show_axis=false)
ges = LScene(fig[2,2], show_axis=false)

Label(fig[1, 1, Top()], "A and B", padding = (0, 0, 1, 0))
Label(fig[1, 2, Top()], "A + B", padding = (0, 0, 1, 0))
Label(fig[2, 1, Top()], "A - B", padding = (0, 0, 1, 0))
Label(fig[2, 2, Top()], "B - A", padding = (0, 0, 1, 0))

render!(ax, s1; color=:red, transparency=true, alpha=0.6)
render!(ax, s2; color=:green)

u = s1+s2
render!(res, u; color=:red)

m = s1-s2
render!(tes, m; color=:red, transparency=true, alpha=0.6)

k = s2-s1
render!(ges, k; color=:red)

display(fig)
set_view(ax, cview)
set_view(res, cview)
set_view(tes, cview)
set_view(ges, cview)

save("composite_sdf.png", fig; px_per_unit=4, update = false)