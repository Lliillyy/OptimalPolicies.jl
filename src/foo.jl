# [src/foo.jl]
"""
    foo(x, y)

Creates a 2-element static array from the scalars `x` and `y`.
"""
function foo(x::Number, y::Number)
    SA[x, y]
end

"""
    U(gam, x)

The potential function for some given potential barrier height `gam`
"""
function U(gam::Number, x::Number)
    gam * (x * x - 1) * (x * x - 1)
end

"""
    curried(gam)

Returns a potential function for the given potential barrier height `gam`
"""
function curried(gam::Number)
    x -> U(gam, x)
end
