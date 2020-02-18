using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using SparseArrays
using DifferentialEquations
#using DiffEqOperators
using Statistics

include("randutils.jl")

struct Metapopulation{T<:AbstractSimpleGraph}
    g::T
    V::Vector{<:Integer} # node capacities
    β::Real # infection rate
    μ::Real # migration rate
end

function normalized_laplacian(g::SimpleGraph)
    L = float(laplacian_matrix(g))
    for i in vertices(g)
        if L[i,i] != 0
            L[i,:] /= L[i,i]
        end
    end
    return L
end

function normalized_laplacian(g::SimpleDiGraph)
    L = float(laplacian_matrix(g))
    for i in vertices(g)
        if L[i,i] != 0
            L[i,:] /= L[i,i]
        end
    end
    return L
end

function metapopulation_gillespie(mp::Metapopulation{SimpleGraph{T}}, s_0, i_0; tmax=100.0, nmax=1000, normalize_migration=true) where T <: Integer
    t = 0.0
    n = 1
    β = mp.β
    μ = mp.μ
    V = mp.V
    N = length(vertices(mp.g))
    ts = Vector{typeof(t)}(undef, nmax)
    ts[1] = t
    ss = Array{typeof(s_0[1]),2}(undef, nmax, length(s_0))
    is = Array{typeof(i_0[1]),2}(undef, nmax, length(i_0))
    ss[1,:] = s_0
    is[1,:] = i_0
    s = copy(s_0)
    i = copy(i_0)
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

function metapopulation_ode(mp::Metapopulation{<:AbstractSimpleGraph{<:Integer}}, s, i; tmax = 100.0)
    L = normalized_laplacian(mp.g)
    # define a special operator for solving the linear part of the equation exactly
    N = nv(mp.g)
    f! = function(dx,x,p,t)
        n = p
        tmp = mp.β*x[1:n].*x[n+1:2*n]./(mp.V).^2
        dx[1:n] = -tmp - μ*L*x[1:n]./mp.V
        dx[n+1:2*n] = tmp - μ*L*x[n+1:2*n]./mp.V
    end
    prob = ODEProblem(f!, vcat(float(s),float(i)), (0.0,tmax), N)
    solve(prob, Tsit5()) # order 5 solver
    # tried to use an exact solver for the linear part, but couldn't get it to work :(
end

function metapopulation_montecarlo(mp::Metapopulation{SimpleGraph{T}}, s_0, i_0; tmax = 100.0, nmax = 1000, nsims = 100, nbins = 100) where T <: Integer
    n = nv(mp.g)
    ts = LinRange(0.0, tmax, nbins)
    s_sum = zeros(nbins, n)
    #s_sumofsquares = zeros(nbins, n)
    i_sum = zeros(nbins, n)
    #i_sumofsquares = zeros(nbins, n)
    for j in 1:nsims
        t, s, i = metapopulation_gillespie(mp, s_0, i_0, tmax=tmax, nmax=nmax)
        l = 1
        for k in 1:nbins
            # Y U NO WORK?
            while l < length(t) && t[l+1] < ts[k]
                l += 1
            end
            s_sum[k,:] += s[l,:]
            i_sum[k,:] += i[l,:]
        end
    end
    mean_s = s_sum/nsims
    mean_i = i_sum/nsims
    #var_s = s_sumofsquares/nsims - mean_s.^2
    #var_i = i_sumofsquares/nsims - mean_i.^2
    return ts, mean_s, mean_i#, sqrt.(var_s), sqrt.(var_i)
end
