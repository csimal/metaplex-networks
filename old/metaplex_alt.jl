using LightGraphs
using LightGraphs.SimpleGraphs
using LinearAlgebra
using DifferentialEquations
using ProgressMeter

include("utils\\randutils.jl")
include("utils\\categorical_tree.jl")
include("utils\\utils.jl")

struct HeterogeneousMetaplex{T1<:AbstractSimpleGraph,T2<:AbstractSimpleGraph}
    g::Vector{T1} # relation networks between individuals
    h::T2 # relation network of components
end


function metaplex_gillespie(mp::HeterogeneousMetaplex, Xi0::BitVector, Xμ0::Vector{Int}, β, D::Vector; nmax=length(Xi0), tmax=100.0, sampling_method=:tree)
    N = nv(mp.g[1]) # number of individuals
    gs = mp.g
    M = nv(mp.h) # number of populations
    h = mp.h
    #od = map(k-> k>0 ? 1/float(k) : 0, outdegree(h))
    #od = outdegree(h)
    od = outdegree(h) .> 0
    Xi = copy(Xi0) # epidemic state of each individual
    Xμ = copy(Xμ0) # location of each individual
    popcounts = zeros(Int,2,M) # number of of individuals susceptible and infected for each component
    for k in 1:N
        if !Xi[k]
            #push!(pops_s[Xμ[k]],k)
            popcounts[1,Xμ[k]] += 1
        else
            #push!(pops_i[Xμ[k]],k)
            popcounts[2,Xμ[k]] += 1
        end
    end
    a = zeros(N) # aggregate reactions for each individual becoming infected and or migrating
    for k in 1:N
        if !Xi[k]
            l = length(filter(i->Xi[i] && Xμ[i]==Xμ[k], inneighbors(gs[Xμ[k]],k)))
            a[k] = β*l # probability of k recieving infection is β time its number of infected incoming neighbors in the same population
            a[k] += D[1]
        else
            a[k] = D[2]
        end
    end

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
        k = rand_categorical(a, a0) # pick a random node
        μ = Xμ[k]
        if !Xi[k]
            if rand()*a[k] < D[1] # node k moves to another cluster
                ν = rand(outneighbors(h,μ))
                Xμ[k] = ν
                popcounts[1,μ] -= 1 # length(pops_s[μ])
                popcounts[1,ν] += 1 #length(pops_s[ν])
                l = length(filter(i->Xi[i] && Xμ[i]==ν, inneighbors(gs[Xμ[k]],k)))
                a[k] = β*l + D[1]
            else # k gets infected
                Xi[k] = true
                popcounts[1,μ] -= 1 #length(pops_s[μ])
                popcounts[2,μ] += 1 #length(pops_i[μ])
                a[k] = D[2]
                for i in Iterators.filter(j->!Xi[j] && Xμ[j]==μ, outneighbors(gs[Xμ[k]],k))
                    a[i] += β
                    #a[i] = D[1] + β*length(filter(i->Xi[i] && Xμ[i]==μ, inneighbors(g,i)))
                end
            end
        else
            ν = rand(outneighbors(h,μ))
            Xμ[k] = ν
            popcounts[2,μ] -= 1
            popcounts[2,ν] += 1 
            for i in Iterators.filter(l->!Xi[l], outneighbors(gs[Xμ[k]],k))
                if Xμ[i] == μ
                    a[i] -= β
                    #a[i] = D[1] + β*length(filter(i->Xi[i] && Xμ[i]==μ, inneighbors(g,i)))
                elseif Xμ[i] == ν
                    a[i] += β
                    #a[i] = D[1] + β*length(filter(i->Xi[i] && Xμ[i]==ν, inneighbors(g,i)))
                end
            end
        end
        a0 = sum(a)
        t += τ
        ts[n+1] = t
        pcs[n+1,:,:] .= popcounts
        n += 1
    end
    return ts[1:n], pcs[1:n,:,:]
end

function metaplex_ode(mp::HeterogeneousMetaplex, Xi0::BitVector, Xμ0::Vector{Int}, β, D::Vector; tmax=100.0, saveat=[])
    N = nv(mp.g[1])
    M = nv(mp.h)
    A = adjacency_matrix.(mp.g)
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
            du[1,i,:] .= .-D[1]*L'*u[1,i,:]
            du[2,i,:] .= .-D[2]*L'*u[2,i,:]
        end
        for μ in 1:M
            infection = β*(u[1,:,μ].*(A[μ]*u[2,:,μ]) )# .+ u[1,:,μ].*u[2,:,μ])
            du[1,:,μ] .-= infection
            du[2,:,μ] .+= infection
        end
    end
    prob = ODEProblem(f!, u0, (0.0,tmax), N)
    sol = solve(prob, Tsit5(), saveat=saveat, abstol=1e-09)
    s = zeros(length(sol.t), N, M)
    i = zeros(length(sol.t), N, M)
    for t in 1:length(sol.t)
        s[t,:,:] = sol.u[t][1,:,:]
        i[t,:,:] = sol.u[t][2,:,:]
    end
    sμ = reshape(sum(s, dims=2), length(sol.t), M)
    iμ = reshape(sum(i, dims=2), length(sol.t), M)
    return sol.t, s, i, sμ, iμ, sol
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
    return ts, Xs_sum/nsims, Xi_sum/nsims, sqrt.(abs.(Ms_sum/(nsims-1))), sqrt.(abs.(Mi_sum/(nsims-1)))
end
