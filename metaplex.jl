using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using DifferentialEquations
using ProgressMeter

include("randutils.jl")
include("categorical_tree.jl")
include("utils.jl")

struct Metaplex{T1<:AbstractSimpleGraph,T2<:AbstractSimpleGraph}
    g::T1 # global relation network between individuals
    h::T2 # relation network of components
end

abstract type EpidemicEvent end
struct InfectionEvent <: EpidemicEvent
    infected::Int
    infector::Int
end

struct MigrationEvent <: EpidemicEvent
    individual::Int
    from::Int
    to::Int
end

function metaplex_gillespie(mp::Metaplex{T1,T2}, Xi0::BitVector, Xμ0::Vector{Int}, β, D::Vector; nmax=length(Xi0), tmax=100.0, sampling_method=:tree) where {T1,T2}
    N = nv(mp.g) # number of individuals
    g = mp.g
    M = nv(mp.h) # number of populations
    h = mp.h
    Xi = copy(Xi0) # epidemic state of each individual
    Xμ = copy(Xμ0) # location of each individual
    popcounts = zeros(Int,2,M) # number of of individuals susceptible and infected for each component
    pops_s = [Set{Int}() for i in 1:M] # Population to individuals without traversing the whole list of individuals
    pops_i = [Set{Int}() for i in 1:M]
    for k in 1:N
        if !Xi[k]
            push!(pops_s[Xμ[k]],k)
            #popcounts[1,Xμ[k]] += 1
        else
            push!(pops_i[Xμ[k]],k)
            #popcounts[2,Xμ[k]] += 1
        end
    end
    popcounts[1,:] = map(length, pops_s)
    popcounts[2,:] = map(length, pops_i)
    a = zeros(N+2*M) # aggregate reactions for each individual becoming infected and individuals migrating from each population
    for k in 1:N
        if !Xi[k]
            l = length(filter(i->Xi[i] && Xμ[i]==Xμ[k], inneighbors(g,k)))
            a[k] = β*l # probability of k recieving infection is β time its number of infected incoming neighbors in the same population
            #nreactions += l
        end
    end
    a[N+1:N+M] .= D[1].*popcounts[1,:].*(outdegree(h) .> 0)
    a[N+M+1:N+2*M] .= D[2].*popcounts[2,:].*(outdegree(h) .> 0)
    if sampling_method == :tree
        a = CategoricalTree(a)
    end
    a0 = sum(a)
    ts = zeros(nmax)
    pcs = zeros(Int,nmax,2,M)
    pcs[1,:,:] .= popcounts
    t = 0.0
    n = 1
    while n < nmax && t < tmax && a0 > 0
        τ = log(1/rand())/a0
        k = rand_categorical(a, a0)
        if k in 1:N # node k becomes infected
            Xi[k] = true
            delete!(pops_s[Xμ[k]], k)
            push!(pops_i[Xμ[k]], k)
            popcounts[1,Xμ[k]] = length(pops_s[Xμ[k]])
            popcounts[2,Xμ[k]] = length(pops_i[Xμ[k]])
            if length(outneighbors(h,Xμ[k])) > 0
                a[N+Xμ[k]] = D[1]*popcounts[1,Xμ[k]]
                a[N+M+Xμ[k]] += D[2]*popcounts[2,Xμ[k]]
            end
            a[k] = 0.0
            for i in Iterators.filter(j->!Xi[j] && Xμ[j]==Xμ[k], outneighbors(g,k))
                a[i] += β
            end
        elseif k in N+1:N+M # a random susceptible node from population k-N moves to another population
            μ = k-N
            ν = rand(outneighbors(h,μ))
            i = rand(pops_s[μ])
            Xμ[i] = ν
            delete!(pops_s[μ], i)
            push!(pops_s[ν], i)
            popcounts[1,μ] = length(pops_s[μ])
            popcounts[1,ν] = length(pops_s[ν])
            a[N+μ] = D[1]*popcounts[1,μ]
            a[N+ν] = D[1]*popcounts[1,ν]*(outdegree(h,ν)>0)
            l = length(filter(j->Xi[j] && Xμ[j]==ν, inneighbors(g,i)))
            a[i] = β*l
        else # a random infected node from k-N-M moves to another population
            μ = k-N-M
            ν = rand(outneighbors(h,μ))
            i = rand(pops_i[μ])
            Xμ[i] = ν
            delete!(pops_i[μ], i)
            push!(pops_i[ν], i)
            for j in Iterators.filter(l->!Xi[l], outneighbors(g,i))
                if Xμ[j] == μ
                    a[j] -= β
                elseif Xμ[j] == ν
                    a[j] += β
                end
            end
            popcounts[2,μ] = length(pops_i[μ])
            popcounts[2,ν] = length(pops_i[ν])
            a[N+M+μ] = D[2]*popcounts[2,μ]
            a[N+M+ν] = D[2]*popcounts[2,ν]
        end
        a0 = sum(a)
        t += τ
        ts[n+1] = t
        pcs[n+1,:,:] .= popcounts
        n += 1
    end
    return ts[1:n], pcs[1:n,:,:]
end

function metaplex_ode(mp::Metaplex, Xi0::BitVector, Xμ0::Vector{Int}, β, D::Vector; tmax=100.0)
    N = nv(mp.g)
    M = nv(mp.h)
    A = adjacency_matrix(mp.g)
    L = normalized_laplacian(mp.h)
    u0 = zeros(2,N,M)
    for i in 1:N
        if !Xi0[i]
            u0[1,i,Xμ0[i]] = 1.0
        else
            u0[2,i,Xμ0[i]] = 1.0
        end
    end
    f! = function(du,u,p,t)
        for i in 1:N
            du[1,i,:] .= .-D[1]*L*u[1,i,:]
            du[2,i,:] .= .-D[2]*L*u[2,i,:]
        end
        for μ in 1:M
            infection = β*u[1,:,μ].*(A*u[2,:,μ])
            du[1,:,μ] .+= .-infection
            du[2,:,μ] .+= infection
        end
    end
    prob = ODEProblem(f!, u0, (0.0,tmax), N)
    solve(prob, Tsit5())
end

function metaplex_montecarlo(mp::Metaplex, Xi0::BitVector, Xμ0::Vector{Int}, β, D::Vector; nmax=length(Xi0), tmax=100.0, nsims=1000, nbins=400)
    ts = LinRange(0.0, tmax, nbins)
    M = nv(mp.h)
    Xs_sum = zeros(nbins, M)
    Ms_sum = zeros(nbins, M)
    Xi_sum = zeros(nbins, M)
    Mi_sum = zeros(nbins, M) # running SSE
    Xn = Vector(undef, M)
    Δn = Vector(undef, M)
    p = Progress(nsims, dt=1.0)
    for n in 1:nsims
        t, pcs = metaplex_gillespie(mp, Xi0, Xμ0, β, D, nmax=nmax, tmax=tmax)
        l = 1
        for k in 1:nbins
            while l < length(t) && t[l+1] < ts[k]
                l += 1
            end
            Xn .= pcs[l,1,:]
            Δn .= Xn .- Xs_sum[k]/(n>1 ? n-1 : 1)
            Xs_sum[k,:] .+= Xn
            Ms_sum[k,:] .+= Δn.*(Xn.-Xs_sum[k]/n)
            Xn .= pcs[l,2,:]
            Δn .= Xn .- Xi_sum[k]/(n>1 ? n-1 : 1)
            Xi_sum[k,:] .+= Xn
            Mi_sum[k,:] .+= Δn.*(Xn.-Xi_sum[k]/n)
        end
        ProgressMeter.next!(p; showvalues= [(:n,n)])
    end
    return ts, Xs_sum/nsims, Xi_sum/nsims, Ms_sum/(nsims-1), Mi_sum/(nsims-1)
end
