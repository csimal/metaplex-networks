using LightGraphs
using LightGraphs.SimpleGraphs

include("metaplex.jl")

N = 100
M = 10

g = random_configuration_model(N*M, fill(6,N*M))
g = complete_graph(N*M)
h = random_configuration_model(M,fill(3,M))
h = path_graph(M)
h = complete_graph(M)
using GraphPlot
gplot(h)

mp = Metaplex(g,h)

Xi0 = falses(N*M)
Xi0[1:100] .= true

Xμ0 = [div(i-1,N)+1 for i in 1:N*M]
Xμ0 = [rand(1:M) for i in 1:N*M]

β = 0.3
D = [0.1, 0.1]
tmax = 3.0

ts, pcs = metaplex_gillespie(mp, Xi0, Xμ0, β, D, nmax = 500000, tmax=tmax)

using Plots
gr()
pyplot()

plot(ts, pcs[:,2,:], label="", line=:steppre)

plot(ts, sum(pcs[:,2,:], dims=2)/(N*M),
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )

plot(ts, vec(sum(pcs[:,:,:], dims=[2,3])), label="")

sol = metaplex_ode(mp, Xi0, Xμ0, β, D, tmax = tmax)

u_s = zeros(length(sol.t), N*M, M)
u_i = zeros(length(sol.t), N*M, M)
for t in 1:length(sol.t)
    u_s[t,:,:] = sol.u[t][1,:,:]
    u_i[t,:,:] = sol.u[t][2,:,:]
end
ui_sum = reshape(sum(u_i, dims=2), length(sol.t), M)

plot(sol.t, reshape(sum(u_i, dims=3), length(sol.t), N*M), label="")

plot(sol.t, ui_sum, label="")

sum(ui_sum, dims=2)/(N*M)

plot(ts, sum(pcs[:,2,:], dims=2)/(N*M),
    label="Gillespie",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals")
plot!(sol.t, sum(ui_sum, dims=2)/(N*M), label="Mean Field")

using ProgressMeter

ts, mean_s, mean_i, var_s, var_i = metaplex_montecarlo(mp, Xi0, Xμ0, β, D, tmax=tmax, nmax=500000)

plot(ts, mean_i, label="")
plot!(sol.t, ui_sum, label="", linestyle=:dash)

plot(ts, sum(mean_i, dims=2)/(N*M), label="Average", legend=:bottomright)
plot!(sol.t, sum(ui_sum, dims=2)/(N*M), label="Mean Field")

sum(mean_i, dims=2)/(N*M)
plot(ts[1:20], sum(mean_i[1:20,:], dims=2)/(N*M), label="Average", legend=:topleft)
plot!(sol.t[1:50], sum(ui_sum, dims=2)[1:50]/(N*M), label="Mean Field")
