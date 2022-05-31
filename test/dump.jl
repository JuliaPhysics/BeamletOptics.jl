"""
    orthogonal(input::Ray)

Returns the normed counter-clockwise orthogonal vector for the input ray direction.
Only works for 2x1 vectors. 
"""
function orthogonal(input::Ray)
    return [-input.dir[2], input.dir[1]] ./ norm(input.dir)
end

"""
    rotate(input::Ray)

Returns the counter-clockwise rotated direction of the input ray.
The vector is rotated θ radians. Length correction is not necessery since rotation maintains length.
Only works for 2x1 vectors. 
"""
function rotate(input::Ray, θ)
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    return R*input.dir
end