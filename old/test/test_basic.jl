using NetworkEpidemics
using Graphs
using Statistics
using Plots
using ColorSchemes
using Random

Random.seed!(2022)

N = 1000
M = 10

h = complete_graph(M)
g = erdos_renyi(N, 0.3)

β = 0.1
D = 0.1

mp = Metapopulation(h, D, SI(β))
mpx = Metaplex(g, h, D, SI(β))

begin
    x0_i = fill(1, N)
    x0_i[1:10] .= 2
    x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]
end

begin
    x0_mp = Array{Int,2}(undef, M, 2)
    x0_mp[2:M,:] .= [100 0]
    x0_mp[1, :] .= [90, 10]
end

begin
    # infected seed randomly dispersed over meta-nodes
    x0_i = fill(1, N)
    x0_i[1:10] .= 2
    x0_μ = rand(1:M, N)
    x0_mp = zeros(Int, M, 2)
    for i in 1:N
        x0_mp[x0_μ[i], x0_i[i]] += 1
    end
end

begin
    tmax = 30.0
    nmax = 100000
    nsims = 1000
    nbins = 200
end

ts_av_mpx, u_av_mpx = average(mpx, [x0_i, x0_μ], tmax=tmax, nmax=nmax, nbins=nbins, nsims=200, progressbar=false)

begin
    plot(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N,
    label="Metaplex",
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    linestyle=:dash
    )
end