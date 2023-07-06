# [src/foo.jl]
"""
    foo(x, y)

Creates a 2-element static array from the scalars `x` and `y`.
"""
function foo(x::Number, y::Number)
    SA[x, y]
end

"""
    U(gam::Number, x::Number)

The potential function for some given potential barrier height `gam`
"""
function U(gam::Number, x::Number)
    gam * (x * x - 1) * (x * x - 1)
end

"""
    curried(gam::Number)

Returns a potential function for the given potential barrier height `gam`
"""
function curried(gam::Number)
    x -> U(gam, x)
end


"""
    chain(target::Function, tune=0.1, init=1, iters=1e5)

The Metropolis algorithm targeting a particular potential function `target`
"""
function chain(target::Function, tune=0.1, init=1, iters=1e5)
    x = init
    xvec = zeros(iters)
    for i in 1:iters
        can = x + rand(Normal(0, tune), 1)
        logA = target(x) - target(can)
        if log(rand(1)) < logA
            x = can
        end
        xvec[i] = x
    end
    return xvec
end