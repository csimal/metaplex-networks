using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using SparseArrays
using DifferentialEquations
#using DiffEqOperators
using Statistics

include("randutils.jl")
include("utils.jl")

struct Metapopulation{T<:AbstractSimpleGraph}
    g::T
    V::Vector{<:Integer} # node capacities
    β::Float64 # infection rate
    μ::Float64 # migration rate
end


"""
    metapopulation_gillespie(mp, s_0, i_0; kwargs)

Perform a single simulation of SI spreading on the metapopulation described by `mp`.

`s0` and `i0` contain respectively the initial numbers of susceptible and infected individuals per node in the metapopulation.

Keyword arguments
  * `nmax=1000`: maximum number of iterations of the algorithm
  * `tmax=100.0`: maximum time of the simulation
  * `normalize_migration=true`: whether or not migration rate should be the same for all nodes (currently, does nothing if set to `false`)
"""
function metapopulation_gillespie(mp::Metapopulation{SimpleGraph{T}}, s0, i0; nmax=1000, tmax=100.0, normalize_migration=true) where T
    t = 0.0
    n = 1
    β = mp.β
    μ = mp.μ
    V = mp.V
    N = length(vertices(mp.g))
    ts = Vector{typeof(t)}(undef, nmax)
    ts[1] = t
    ss = Array{typeof(s0[1]),2}(undef, nmax, length(s0))
    is = Array{typeof(i0[1]),2}(undef, nmax, length(i0))
    ss[1,:] = s0
    is[1,:] = i0
    s = copy(s0)
    i = copy(i0)
    d = degree(mp.g) .> 0 # which nodes have neighbors
    a = Vector{Float64}(undef, 3*N)
    a[1:N] = β*s.*i./(V.^2)
    a[N+1:2*N] = μ*s./V.*d # NB. aggregated probabilities of migration FROM each node
    a[2*N+1:3*N] = μ*i./V.*d # isolated nodes have migration probability of 0

    a0 = sum(a)
    while t < tmax && n < nmax && a0 > 0
        #println(a)
        τ = log(1/rand())/a0
        k = rand_categorical(a)
        if k <= N # infection in node k
            s[k] -= 1
            i[k] += 1
            a[k] = β*s[k]*i[k]/(V[k]^2)
            if d[k]
                a[N+k] = μ*s[k]/V[k]
                a[2*N+k] = μ*i[k]/V[k]
            end
        elseif k <= 2*N # migration of susceptible individual
            u = k-N
            v = rand(neighbors(mp.g,u)) # NB. will fail if there are no neighbors
            s[u] -= 1
            s[v] += 1
            a[u] = β*s[u]*i[u]/(V[u]^2)
            a[N+u] = μ*s[u]/V[u]
            a[v] = β*s[v]*i[v]/(V[v]^2)
            a[N+v] = μ*s[v]/V[v]
        else # migration of infected individual
            u = k-2*N
            v = rand(neighbors(mp.g,u)) # NB. will fail if there are no neighbors
            i[u] -= 1
            i[v] += 1
            a[u] = β*s[u]*i[u]/(V[u]^2)
            a[2*N+u] = μ*i[u]/V[u]
            a[v] = β*s[v]*i[v]/(V[v]^2)
            a[2*N+v] = μ*i[v]/V[v]
        end
        a0 = sum(a)
        t += τ
        ts[n+1] = t
        ss[n+1,:] = s
        is[n+1,:] = i
        n += 1
    end
    return ts[1:n], ss[1:n,:], is[1:n,:]
end

"""
    metapopulation_ode(mp, s0, i0; kwargs)

Integrate the Mean Field equations for SI spreading on the metapopulation described by `mp`.

`s0` and `i0` contain respectively the initial numbers of susceptible and infected individuals per node.

Keyword arguments:
  * `tmax=100.0`: the maximum time of the simulation
"""
function metapopulation_ode(mp::Metapopulation{<:AbstractSimpleGraph{<:Integer}}, s0, i0; tmax = 100.0)
    L = normalized_laplacian(mp.g)
    # define a special operator for solving the linear part of the equation exactly
    N = nv(mp.g)
    f! = function(dx,x,p,t)
        n = p
        tmp = mp.β*x[1:n].*x[n+1:2*n]./(mp.V).^2
        dx[1:n] = -tmp - μ*L'*x[1:n]./mp.V
        dx[n+1:2*n] = tmp - μ*L'*x[n+1:2*n]./mp.V
    end
    prob = ODEProblem(f!, vcat(float(s0),float(i0)), (0.0,tmax), N)
    solve(prob, Tsit5()) # order 5 solver
    # tried to use an exact solver for the linear part, but couldn't get it to work :(
end

"""
    metapopulation_montecarlo(mp, s0, i0; kwargs)

Perform multiple simulations of SI spreading on the metapopulation described by `mp`, returning the average curve.

`s0` and `i0` contain respectively the initial numbers of susceptible and infected individuals in each node.

Keyword arguments:
  * `nmax=1000`: maximum number of iterations of each simulation
  * `tmax=100.0`: maximum time of the simulations
  * `nsims=100`: number of simulations to perform
  * `nbins`: the number of time steps at which the average is computed
"""
function metapopulation_montecarlo(mp::Metapopulation{SimpleGraph{T}}, s0, i0; nmax = 1000, tmax = 100.0, nsims = 100, nbins = 100) where T <: Integer
    N = nv(mp.g)
    ts = LinRange(0.0, tmax, nbins)
    Xs_sum = zeros(nbins, N)
    Ms_sum = zeros(nbins, N)
    Xi_sum = zeros(nbins, N)
    Mi_sum = zeros(nbins, N)
    Xn = Vector(undef, N)
    Δn = Vector(undef, N)
    for n in 1:nsims
        t, s, i = metapopulation_gillespie(mp, s0, i0, tmax=tmax, nmax=nmax)
        l = 1
        for k in 1:nbins
            while l < length(t) && t[l+1] < ts[k]
                l += 1
            end
            Xn .= s[l,:]
            Δn .= Xn .- Xs_sum[k,:]/(n>1 ? n-1 : 1)
            Xs_sum[k,:] .+= Xn
            Ms_sum[k] .+= Δn.*(Xn.-Xs_sum[k]/n)
            Xn .= i[l,:]
            Δn .= Xn .- Xi_sum[k,:]/(n>1 ? n-1 : 1)
            Xi_sum[k,:] .+= Xn
            Mi_sum[k,:] .+= Δn.*(Xn.-Xi_sum[k]/n)
        end
    end
    return ts, Xs_sum/n, Xi_sum/n, Ms_sum/(n-1), Mi_sum/(n-1)
end
