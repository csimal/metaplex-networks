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

function metapopulation_gillespie(mp::Metapopulation{SimpleGraph{T}}, s_0, i_0; tmax=100.0, nmax=1000) where T <: Integer
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
        ts[n+1] = t
        ss[n+1,:] = s
        is[n+1,:] = i
        t += τ
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

function metapopulation_gillespie_montecarlo(mp::Metapopulation{SimpleGraph{<:Integer}}, s_0, i_0; tmax = 100.0, nmax = 1000, nsims = 100, nbins = 100)
    n = length(vertices(mp.g))
    ts = LinRange(0.0, tmax, nbins)
    s_sum = zeros(nbins, n)
    #s_sumofsquares = zeros(nbins, n)
    i_sum = zeros(nbins, n)
    #i_sumofsquares = zeros(nbins, n)
    for j in 1:nsims
        t, s, i = metapopulation_gillespie(mp, s_0, i_0, tmax=tmax, nmax=nmax)
        l = 1
        for k in 1:nbins
            while ts[k] > t[l] && l < length(t)
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

N = 10
P = 1000
volume = 1
frac_infected = 0.01
ninfected = Int(frac_infected*P)
g = complete_graph(N)
g = random_configuration_model(N,fill(3,N))
V = fill(volume, N)
β = 0.1
μ = 0.3
mp = Metapopulation(g,V,β,μ)
s0 = fill(P, N)
i0 = zeros(Int, N)
s0[1] = P-ninfected
i0[1] = ninfected
Ptot = sum(s0) + sum(i0)
ts, s, i = metapopulation_gillespie(mp,s0,i0, nmax = 1000000, tmax=0.8)

using Plots
plot(ts, sum(i, dims=2)/Ptot, label="", line=:steppre)

plot(ts,i/P, line=:steppre, label="")
plot(ts,s/Ptot, line=:steppre)

sol = metapopulation_ode(mp, s0, i0, tmax=0.8)
plot(sol, vars = collect(N+1:2*N), label="")
plot!(ts,i, line=:steppre, label="")

plot(ts,sum(i+s, dims=2), label="")
all(x->x==1001, sum(i+s, dims=2))

tsmc, mean_s, mean_i = metapopulation_gillespie_montecarlo(mp,s0,i0,tmax=2.0, nmax=5000, nbins=200)

plot(tsmc, mean_i, line=:steppre)
plot!(sol, vars = collect(N+1:2*N))

function transient_time_montecarlo(mp::Metapopulation, s0, i0; tmax=100.0, nmax=1000, nsims = 100)
    totalpop = sum(s0) + sum(i0)
    tts = zeros(nsims)
    for l in 1:nsims
        t, s, i = metapopulation_gillespie(mp, s0, i0, tmax=tmax, nmax=nmax)
        k = 1
        tot = sum(i[k,:])
        while tot < 0.9*totalpop && k < length(t)
            k += 1
            tot += sum(i[k,:])
        end
        tts[l] = t[k]
    end
    return mean(tts)
end

μs = LinRange(0.0,1.0, 200)
ttime = zeros(200)
for k in 1:200
    mp = Metapopulation(g, β, μs[k])
    ttime[k] = transient_time_montecarlo(mp,s0,i0,tmax = 1.0, nmax=5000, nsims=200)
end
scatter(μs, ttime, label="", xlabel="\\mu", ylabel="\\tau")
savefig("transientvsmu.png")

βs = LinRange(0.0,1.0,200)

ttimes = zeros(200,200)

for k in 1:200, l in 1:200
    mp = Metapopulation(g, βs[k], μs[l])
    ttimes[k,l] = transient_time_montecarlo(mp,s0,i0,tmax=1.0,nmax=5000,nsims=200)
end
