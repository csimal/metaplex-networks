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
h = london_transport()
h = erdos_renyi(M, 0.3)

M = nv(h)
N = M*1000
graphplot(h)

begin
    p = sortperm(degree(h))
    eig = eigen(Array(normalized_laplacian(h)))
    V = Array{Float64,2}(undef, nv(h), nv(h))
    for i in 1:nv(h)
        V[:,i] .= eig.vectors[p,i]
    end
    # for each node, find the eigenvector whose entry for that node is highest in magnitude
    eigenmode = Vector{Int}(undef, nv(h))
    for i in 1:nv(h)
        eigenmode[i] = findmax(abs.(eig.vectors[i,:]))[2]
    end
end

histogram(degree(h), label="")
scatter!(degree_histogram(h), label="")

scatter(eig.values, label="")

heatmap(V', color=cgrad(:balance))

bar(eigenmode, label="")

β = 0.39
γ = 0.4
D = 0.1

critical_k = (N / M) * γ / β


ks = fill(50, M)
ks_ = copy(ks)
ks_[1] = 400

gs = [random_regular_graph(N, k) for k in ks]
mpx = HeterogeneousMetaplex(gs, h, D, SIS(β,γ))

mp = HeterogeneousMetapopulation(h, N, ks, D, SIS(β,γ))
mp_ = HeterogeneousMetapopulation(h, N, ks_, D, SIS(β,γ))

x0_i = fill(1, N)
x0_i[1:10] .= 2

x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]
x0_μ = [mod(i,M)+1 for i in 1:N]

x0 = Array{Int,2}(undef, M, 2)
x0[2:M,:] .= [1000 0]
x0[1,:] .= [980, 20]

tmax = 100.0
nmax = 1000000
nsims = 100
nbins = 200

#ts_mpx, u_mpx = gillespie(mpx, [x0_i,x0_μ], nmax=nmax, tmax=tmax)
ts_mp, u_mp = gillespie(mp, x0, nmax=nmax, tmax=tmax)

u_mp

#plot(ts_mpx, u_mpx[2], label="")
plot(ts_mp, u_mp[2], label="")

#ts_av_mpx, u_av_mpx = average(mpx, [x0_i, x0_μ], nmax=nmax, tmax=tmax, nbins=nbins, nsims=nsims)
#ts_mf_mpx, u_mf_mpx = meanfield(mpx, [x0_i,x0_μ], tmax=tmax)

ts_av_mp, u_av_mp = average(mp, x0, nmax=nmax, tmax=tmax, nbins=nbins, nsims=nsims, progressbar=true)
ts_mf_mp, u_mf_mp = meanfield(mp, x0, tmax=200.0)
ts_mf_mp_, u_mf_mp_ = meanfield(mp_, x0, tmax=200.0)
#plot(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N, label="Average", ylims=(0.0,1.0))
#plot!(ts_mf_mpx, sum(u_mf_mpx[2], dims=2)/N,  label="Meanfield")

plot(ts_av_mp, sum(u_av_mp[2], dims=2)/N, label="Average")
plot!(ts_mf_mp, sum(u_mf_mp[2], dims=2)/N, label="Meanfield")

plot(ts_mf_mp, sum(u_mf_mp[2], dims=2)/N, label="Meanfield", legend=:topleft)
plot!(ts_mf_mp_, sum(u_mf_mp_[2], dims=2)/N, label="Meanfield", legend=:topleft)

plot(ts_av_mp, u_av_mp[2], label="")
plot(ts_mf_mp, u_mf_mp[2], label="")

#plot(ts_av_mpx, u_av_mpx[2], label="")
#plot(ts_mf_mpx, u_mf_mpx[2], label="")

v = inv(eig.vectors)
ψ_mf_mp = copy(u_mf_mp[2])
for i in 1:size(u_mf_mp)[1]
   ψ_mf_mp[i,:] .= v * u_mf_mp[2][i,:]
end
plot(ts_mf_mp, ψ_mf_mp, label="")

plot(ts_mf_mp, (eig.vectors * u_mf_mp[2]')', label="")