using NetworkEpidemics
using LightGraphs

N = 10
g = complete_graph(N)

β = 0.3
D = 0.1

mp = Metapopulation(g, D, SI(β))

x0 = Array{Int,2}(undef, N,2)
x0[2:N,:] .= [100 0]
x0[1, :] .= [90, 10]

ts, output = gillespie(mp, x0, nmax=2000)

using Plots

plot(ts, output[2], line=:steppre, label="")
