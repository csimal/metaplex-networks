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

function SI(α::Real)
    rs = ReactionSystem(2, 1, [α],
    [X -> X[1]*X[2]],
    [[(1,-1),(2,1)]],
    [[1]]
    )
end

function gillespie_algorithm(rs::ReactionSystem, Xi, tmax = 100.0, nmax = 1000)
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

si = SI(0.5)
si.h[1]([10,10])

ts, Xs = gillespie_algorithm(si, [100, 1])

using Plots

plot(ts, Xs[:,1])
scatter!(ts, Xs[:,1], markersize=1)
plot!(ts, [Xs[i][2] for i in 1:length(Xs)])
