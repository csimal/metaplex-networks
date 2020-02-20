using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using SparseArrays
using DifferentialEquations

include("randutils.jl")
include("categorical_tree.jl")

"""
    contact_process_gillespie(g, Xi, β; kwargs)

Simulate SI spreading on network `g` using the Gillespie algorithm.

`g` is an `AbstractSimpleGraph{<:Integer}`, `β` is the infection rate and `Xi` is a bit vector whose true entries denote initially infected nodes. Additional arguments are
  * `nmax=1000`: the maximum number of iterations of the algorithm
  * `tmax=100.0`: the maximum time allowed
  * `sampling_method=:tree`: The algorithm for sampling which node gets infected at a given time step. Choose from `:array`, `:sparse` or `:tree`. Default is `:tree` which scales better with the number of nodes (log(N) vs. N^2).
  * `record=false`: Whether or not to record which node gets infected and by whom. If `true`, return a vector of pairs `(k,j)` describing which node `k` got infected and `j` its neighbor that infected it.
"""
function contact_process_gillespie(g::AbstractSimpleGraph{<:Integer}, Xi::BitVector, β::Real; nmax=length(Xi), tmax=100, sampling_method=:tree, record=false)
    N = nv(g)
    X = copy(Xi)
    #nreactions::Int64 = 0
    if sampling_method == :tree
        a = zeros(N)
    elseif sampling_method == :sparse
        a = spzeros(N)
    else
        a = zeros(N)
    end
    for k in 1:N
        if !X[k]
            l = length(filter(i->X[i], inneighbors(g,k)))
            a[k] = β*l# probability of k recieving infection is β time its number of infected incoming neighbors
            #nreactions += l
        end
    end
    if sampling_method == :tree # use a specialized binary tree
        a = CategoricalTree(a)
    end
    a0 = sum(a)
    #a0 = nreactions*β
    Es = Vector{Tuple{Int64,Int64}}(undef, nmax-1) # which node gets infected, and by whom
    ts = zeros(nmax)
    t = 0.0
    n = 1
    while n < nmax && t < tmax && a0 > 0
        τ = log(1/rand())/a0
        k = rand_categorical(a, a0) # which node gets infected
        X[k] = true
        #nreactions -= Int(floor(a[k]/β))
        a[k] = 0.0
        for i in filter(j->!X[j],outneighbors(g,k))
            a[i] += β
            #nreactions += 1
        end
        if sampling_method == :sparse
            SparseArrays.dropstored!(a,k) # O(1) cost?
        end
        #a0 = β*nreactions
        a0 = sum(a)
        t += τ
        ts[n+1] = t
        j = rand(filter(l->X[l],inneighbors(g,k)))
        Es[n] = (k,j)
        n += 1
    end
    if record
        return ts[1:n], sum(Xi) .+ collect(0:n-1), Es[1:(n-1)]
    else
        return ts[1:n], sum(Xi) .+ collect(0:n-1)
    end
end

"""
    contact_process_ode(g, Xi, β; kwargs)

Integrate the mean field approximation for SI spreading on network `g`.

`Xi` is the initial condition, a vector whose elements denote infected nodes by 1 and susceptible by 0. `β` is the infection rate.

Keyword arguments:
  * `tmax = 100.0`: maximum time of the simulation.
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

Keyword arguments:
  * `nmax=1000`: maximum number of iterations of each simulation
  * `tmax=100.0`: maximum time of the simulations
  * `nsims=100`: number of simulations to perform
  * `nbins`: the number of time steps at which the average is computed
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
            Xn = X[l]
            Δn = Xn - X_sum[k]/(n>1 ? n-1 : 1)
            X_sum[k] += Xn
            M_sum[k] += Δn*(Xn-X_sum[k]/n)
        end
    end
    return ts, X_sum/nsims, M_sum/(nsims-1)
end
