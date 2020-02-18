using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using DifferentialEquations

include("randutils.jl")
include("categorical_tree.jl")

"""
    contact_process_gillespie(g, Xi, β; kwargs)

Simulate SI spreading on network `g` using the Gillespie algorithm.

`g` is an `AbstractSimpleGraph{<:Integer}`, `β` is the infection rate and `Xi` is a bit vector whose true entries denote initially infected nodes. Additional arguments are
`nmax=1000`: the maximum number of iterations of the algorithm
`tmax=100.0`: the maximum time allowed
"""
function contact_process_gillespie(g::SimpleDiGraph{<:Integer}, Xi::BitVector, β::Real; nmax=1000, tmax=100)
    N = nv(g)
    E = collect(edges(g))
    X = copy(Xi)
    a = zeros(ne(g))
    for k in 1:ne(g)
        i = src(E[k])
        j = dst(E[k])
        a[k] = β*X[i]*!X[j] # i infected and j susceptible
    end
    #ct = CategoricalTree{Float64}(a)
    a0 = sum(a)
    #a0 = sum(ct)
    Es = Vector{Int64}(undef, nmax-1) # record which edges transmit infection
    Xs = falses(nmax, N)
    Xs[1,:] .= Xi
    X = copy(Xi)
    ts = zeros(nmax)
    t = 0.0
    n = 1
    while n < nmax && t < tmax && !all(Xi) && a0 > 0
        τ = log(1/rand())/a0
        k = rand_categorical(a)
        #k = rand_categorical(ct)
        j = dst(E[k])
        X[j] = true
        for e in 1:ne(g)
            i = dst(E[e])
            if src(E[e]) == j
                a[e] = β*!X[i]
                #ct[e] = β*!X[i]
            end
            if i == j
                a[e] = 0.0
                #ct[e] = 0.0
            end
        end
        a0 = sum(a)
        #a0 = sum(ct)
        t += τ
        ts[n+1] = t
        Xs[n+1,:] = X
        Es[n] = k
        n += 1
    end
    return ts[1:n], Xs[1:n,:], Es[1:(n-1)]
end

function contact_process_gillespie(g::SimpleGraph{<:Integer}, Xi::BitVector, β::Real; nmax=1000, tmax = 100.0)
    return contact_process_gillespie(SimpleDiGraph(g), Xi, β, nmax=nmax, tmax=tmax) # needed to have edges in both directions
end

"""
    contact_process_ode(g, Xi, β; kwargs)

Integrate the mean field approximation for SI spreading on network `g`.

`Xi` is the initial condition, a vector whose elements denote infected nodes by 1 and susceptible by 0. `β` is the infection rate.

Keyword arguments:
`tmax = 100.0`: maximum time of the simulation.
"""
function contact_process_ode(g::AbstractSimpleGraph{<:Integer}, Xi, β::Real; tmax=100.0)
    A = adjacency_matrix(g)
    u0 = float(Xi)
    f! = function(du,u,p,t)
        du .= β*(1.0 .- u).*(A*u)
    end
    prob = ODEProblem(f!, u0, (0.0,tmax))
    return solve(prob, Tsit5())
end

"""
    contact_process_montecarlo(g, Xi, β; kwargs)

Compute the average and variance of trajectories given by the gillespie algorithm over `nsims` simulations.
"""
function contact_process_montecarlo(g::AbstractSimpleGraph{<:Integer}, Xi::BitVector, β::Real; nmax = 1000, tmax = 100.0, nsims=1000, nbins = 100)
    ts = LinRange(0.0, tmax, nbins)
    X_sum = zeros(nbins)
    M_sum = zeros(nbins) # running SSE
    for n in 1:nsims
        t, X = contact_process_gillespie(g, Xi, β, nmax=nmax, tmax=tmax)
        l = 1
        for k in 1:nbins
            while l < length(t) && t[l+1] < ts[k]
                l += 1
            end
            Xn = sum(X[l,:])
            Δn = Xn - X_sum[k]/(n>1 ? n-1 : 1)
            X_sum[k] += Xn
            M_sum[k] += Δn*(Xn-X_sum[k]/n)
        end
    end
    return ts, X_sum/nsims, M_sum/(nsims-1)
end
