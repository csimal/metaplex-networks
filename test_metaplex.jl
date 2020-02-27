using LightGraphs
using LightGraphs.SimpleGraphs

include("metaplex.jl")

N = 1000
M = 10

g = random_configuration_model(N*M, fill(6,N*M))
h = random_configuration_model(M,fill(3,M))

using GraphPlot
gplot(h)

mp = Metaplex(g,h)

Xi0 = falses(N*M)
Xi0[1:10] .= true

Xμ0 = [div(i-1,N)+1 for i in 1:N*M]
Xμ0 = [rand(1:M) for i in 1:N*M]

β = 0.3
D = [0.1, 0.1]

ts, pcs = metaplex_gillespie(mp, Xi0, Xμ0, β, D, nmax = 500000, tmax=2000.0)
pcs[:,2,:]

using Plots
gr()
pyplot()
plot(ts, pcs[:,2,:], label="", line=:steppre)

plot(ts, pcs[:,1,:], label="", line=:steppre)

plot(ts, sum(pcs[:,2,:], dims=2)/(N*M),
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )

plot(ts, vec(sum(pcs[:,:,:], dims=[2,3])), label="")

pcs[:,1,:]

sol = metaplex_ode(mp, Xi0, Xμ0, β, D, tmax = 500.0)

us = zeros(length(sol.t), 2, M)
for t in 1:length(sol.t)
    us[t,:,:] = sum(sol.u[t], dims=2)
end
us
u_sum = reshape(sum(us[:,:,:], dims=3), length(sol.t), 2)
vec(u_sum[:,1])
plot(sol.t, us[:,2,:], label="")

plot(sol.t, vec(u_sum[:,2])/(N*M), label="")

plot(ts, sum(pcs[:,2,:], dims=2)/(N*M),
    label="Gillespie",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals")
plot!(sol.t, vec(u_sum[:,2])/(N*M), label="Mean Field")

using ProgressMeter

ts, mean_s, mean_i, var_s, var_i = metaplex_montecarlo(mp, Xi0, Xμ0, β, D, tmax=500.0, nmax=500000)

plot(ts, mean_i, label="")
plot!(sol.t, us[:,2,:], label="", linestyle=:dash)

plot(ts, sum(mean_i, dims=2)/(N*M), label="Average", legend=:bottomright)
plot!(sol.t, vec(u_sum[:,2])/(N*M), label="Mean Field")
