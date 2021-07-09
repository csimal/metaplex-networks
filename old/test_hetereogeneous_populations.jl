using Base: permutecols!!
using NetworkEpidemics
using LightGraphs
using Plots
using Random
using Statistics
using LinearAlgebra
using GraphRecipes

cd("metaplex-networks/old")

include("networks/empirical_networks.jl")

Random.seed!(2021);

N = 1000
M = 20

h = contiguous_usa()

graphplot(h)

begin
    p = sortperm(degree(h))
    eig = eigen(Array(laplacian_matrix(h)))
    V = Array{Float64,2}(undef, nv(h), nv(h))
    for i in 1:nv(h)
        V[i,:] .= eig.vectors[i,p]
    end
    
end

histogram(degree(h))
scatter!(degree_histogram(h))

heatmap(V)

h = erdos_renyi(M, 0.3)

ks = fill(10, M)
ks[1] = 40

gs = [random_regular_graph(N, k) for k in ks]

β = 0.3
γ = 0.1

mpx = HeterogeneousMetaplex(gs, h, 0.1, SIS(β,γ))

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = [div(i-1, div(N,M))+1 for i in 1:N]

tmax = 120.0
nmax = 500000
nsims = 500
nbins = 200

ts_mpx, u_mpx = gillespie(mpx, [x0_i,x0_μ], nmax=nmax, tmax=tmax)

plot(ts_mpx, u_mpx[2])

ts_av_mpx, u_av_mpx = average(mpx, [x0_i, x0_μ], nmax=nmax, tmax=tmax, nbins=nbins, nsims=nsims)
ts_mf_mpx, u_mf_mpx = meanfield(mpx, [x0_i,x0_μ], tmax=tmax)

plot(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N, label="Average")
plot!(ts_mf_mpx, sum(u_mf_mpx[2], dims=2)/N,  label="Meanfield")
