using Statistics
using Random

struct ReactionSystem
    N::Integer
    M::Integer
    c::Vector{T} where T <: Real
    h::Vector{Function}
    update_X::Vector{Vector{Tuple{Integer,Integer}}}
    update_a::Vector{Vector{Integer}}
end

function gillespie_algorithm(rs::ReactionSystem, Xi, tmax = 100.0, nmax = 500)
    t = 0.0
    n = 0
    ts = Vector{typeof(t)}(undef, nmax)
    X = copy(Xi)
    Xs = Array{typeof(Xi[1]),2}(undef, nmax, rs.N)
    Xs[1,:] = X
    a = [rs.c[ν]*(rs.h[ν](X)) for ν in 1:rs.M]
    a0 = sum(a)
    while t < tmax && n < nmax && a0 > 0
        τ = log(1/rand())/a0
        r = rand()
        μ = 1
        s = a[μ]
        while r*a0 > s
            μ += 1
            s += a[μ]
        end
        for (i,j) in rs.update_X[μ]
            X[i] += j
        end
        for ν in rs.update_a[μ]
            a[ν] = rs.c[ν]*(rs.h[ν](X))
        end
        a0 = sum(a) # recompute the entire sum to avoid rounding errors
        t += τ
        n += 1
        ts[n+1] = t
        Xs[n+1,:] = X
    end
    return ts[1:n], Xs[1:n,:]
end

function SI(α::Real)
    rs = ReactionSystem(2, 1, [α],
    [X -> X[1]*X[2]],
    [[(1,-1),(2,1)]],
    [[1]]
    )
end

function LotkaVolterra(c)
    rs = ReactionSystem(3,3,c,
    [X -> X[1]*X[3],
     X -> X[1]*X[2],
     X -> X[2]],
    [[(1,1)],[(1,-1),(2,1)],[(2,-1)]],
    [[1],[1,2],[2]]
    )
end

function Brusselator(c)
    rs = ReactionSystem(4,4,c
    [X -> X[3],
     X -> X[4]*X[1],
     X -> X[2]*X[1](X[1]-1)/2,
     X -> X[2]],
     [[(1,1)],[(1,-1),(2,1)],[(1,-2),(2,2)],[(1,-1)]],
     [[2,3,4],[2,3,4],[2,3,4],[2,3,4]]
    )
end

si = SI(0.5)

ts, Xs = gillespie_algorithm(si, [100,1])

using Plots

plot(ts, Xs[:,1], line=:steppre)
scatter!(ts, Xs[:,1], markersize=1)
plot!(ts, 100*exp.(-0.5*ts))
