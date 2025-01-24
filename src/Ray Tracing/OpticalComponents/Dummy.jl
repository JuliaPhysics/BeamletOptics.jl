struct DummyObject{T, S <: AbstractShape{T}} <: AbstractObject{T, S}
    shape::S
end

set_new_origin3d!(d::DummyObject) = set_new_origin3d!(d.shape)
intersect3d(::DummyObject, ::AbstractRay) = nothing
interact3d(::AbstractSystem, ::DummyObject, ::AbstractBeam, ::AbstractRay) = nothing

MeshDummy(loadpath::String) = DummyObject(Mesh(load(loadpath)))

render_object!(axis, dummy::DummyObject) = render_dummy_mesh!(axis, dummy)

render_dummy_mesh!(::Any, ::DummyObject; transparency = false) = nothing
