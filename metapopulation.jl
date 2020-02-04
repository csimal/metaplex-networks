using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using Statistics
using Random

struct Metapopulation
    g::AbstractSimpleGraph
    β::Real # infection rate
    μ::Real # migration rate
end

function normalized_laplacian(g::SimpleGraph)
    L = laplacian_matrix(g)
    k = degree(g)
    for i in vertices(g)
        if k[i] != 0
            L[i,:] /= k[i]
        end
    end
    return L
end

function normalized_laplacian(g::SimpleDiGraph)
    L = laplacian_matrix(g)
    k = outdegree(g)
    for i in vertices(g)
        if k[i] != 0
            L[i,:] /= k[i]
        end
    end
    return L
end

function rand_mbernoulli(x::Vector{T}) where T <: Real
    r = rand()
    s = sum(x)
    i = 1
    t = x[i]
    while r*s > t
        i += 1
        t += x[i]
    end
    return i
end

function gillespie_metapopulation(mp::Metapopulation, s, i; tmax=100.0, nmax=1000)
    t = 0.0
    n = 0
    L = normalized_laplacian(mp.g)
    ts = Vector{typeof(t)}(undef, nmax)

    while t < tmax && n < nmax && a0 > 0
        τ = log(1/rand())/a0


        t += τ
        n += 1
    end
end
