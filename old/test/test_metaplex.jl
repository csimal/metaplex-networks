using Graphs
using Graphs.SimpleGraphs

include("SI\\metaplex.jl")

N = 100
M = 10

g = random_configuration_model(N*M, fill(6,N*M))
g = complete_graph(N*M)
h = random_configuration_model(M,fill(3,M))
h = erdos_renyi(M, 0.5)
h = path_graph(M)
h = complete_graph(M)

using GraphPlot
gplot(h)

mpx = Metaplex(g,h)

Xi0 = falses(N*M)
Xi0[1:10] .= true

Xμ0 = [div(i-1,N)+1 for i in 1:N*M]
Xμ0 = [rand(1:M) for i in 1:N*M]

β = 0.2
D = [0.1, 0.1]
tmax = 3.0

@time t, pcs = metaplex_gillespie(mpx, Xi0, Xμ0, β, D, nmax = 500000, tmax=tmax)

using Plots
gr()
pyplot()

plot(t, pcs[:,2,:], label="", line=:steppre)

plot(t, sum(pcs[:,2,:], dims=2)/(N*M),
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )

plot(t, vec(sum(pcs[:,:,:], dims=[2,3])), label="")

t_mf, s_mf, i_mf, sμ_mf, iμ_mf, sol = metaplex_ode(mpx, Xi0, Xμ0, β, D, tmax = tmax)


plot(t_mf, reshape(sum(i_mf, dims=3), length(sol.t), N*M), label="")

plot(t_mf, iμ_mf, label="")

sum(iμ_mf, dims=2)/(N*M)

plot(t, sum(pcs[:,2,:], dims=2)/(N*M),
    label="Gillespie",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals")
plot!(t_mf, sum(iμ_mf, dims=2)/(N*M), label="Mean Field")

t_mc, s_mc, i_mc, sd_s_mc, sd_i_mc = metaplex_montecarlo(mpx, Xi0, Xμ0, β, D, tmax=tmax, nmax=500000)

plot(t_mc, i_mc, label="")
plot!(t_mf, iμ_mf, label="", linestyle=:dash)

plot(t_mc, sum(i_mc, dims=2)/(N*M), label="Average", legend=:bottomright)
plot!(t_mf, sum(iμ_mf, dims=2)/(N*M), label="Mean Field")

sum(i_mf, dims=2)/(N*M)
plot(t_mc[1:20], sum(i_mc[1:20,:], dims=2)/(N*M), label="Average", legend=:topleft)
plot!(t_mf[1:50], sum(iμ_mf, dims=2)[1:50]/(N*M), label="Mean Field")
