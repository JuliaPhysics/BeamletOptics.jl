"""
    Nullable{T}

An alias which results in `Union{T, Nothing}` to provide a shorter notation for struct
fields which can containing nothing.
"""
const Nullable{T} = Union{T,Nothing} where {T}

"""
    NullableVector{T}

An alias which results in `Union{Vector{T}, Nothing}` to provide a shorter notation for struct
fields which can containing nothing.
"""
const NullableVector{T} = Union{Vector{T},Nothing} where {T}

include("LinearAlgebraUtils.jl")
include("OpticUtils.jl")
include("RefractiveIndexUtils.jl")

"""
    ToDO
"""
function find_zero_bisection(f, a, b; tol=1e-10, max_iter=1000)
    fa = f(a)
    fb = f(b)
    if sign(fa) == sign(fb)
        error("Bisection requires a sign change: f(a)=$(fa), f(b)=$(fb)")
    end
    for i in 1:max_iter
        mid = (a + b) / 2
        fmid = f(mid)
        if abs(fmid) < tol
            return mid
        end
        if sign(fa) == sign(fmid)
            a = mid
            fa = fmid
        else
            b = mid
            fb = fmid
        end
    end
    error("Bisection did not converge after $max_iter iterations")
end

