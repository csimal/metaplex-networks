using LightGraphs
using LightGraphs.SimpleGraphs
using DifferentialEquations
using Plots

include("metapopulation.jl")

N = 100
P = 1000
volume = 1
frac_infected = 0.01
ninfected = Int(frac_infected*P)
g = complete_graph(N)
g = random_configuration_model(N,fill(3,N))
g = path_graph(N)
V = fill(volume, N)
β = 0.0
μ = 0.5
tmax = 10.0

mp = Metapopulation(g,V,β,μ)
s0 = fill(P, N)
i0 = fill(ninfected, N)
i0 = zeros(Int, N)
s0[1] = P-ninfected
i0[1] = ninfected
Ptot = sum(s0) + sum(i0)
ts, s, i = metapopulation_gillespie(mp,s0,i0, nmax = 1000000, tmax=tmax)


plot(ts, sum(i, dims=2)/Ptot,
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )

plot(ts,i/P,
    line=:steppre,
    label="",
    xlabel = "Time",
    ylabel = "Normalized #infected individuals"
    )
plot(ts,s/Ptot, line=:steppre)

sol = metapopulation_ode(mp, s0, i0, tmax=tmax)

plot(sol, vars = collect(N+1:2*N),
    label="",
    xlabel="Time",
    ylabel="#infected individuals"
    )
plot!(ts,i, line=:steppre, label="")

us = hcat([sol.u[t][N+1:2*N] for t in 1:length(sol.t)]...)
sum(us, dims=1)
plot(sol.t, sum(us, dims=1)')

tsmc, mean_s, mean_i = metapopulation_montecarlo(mp,s0,i0,tmax=tmax, nmax=150000, nbins=200)

plot(tsmc, mean_i,
    label="",
    xlabel="Time",
    ylabel="#infected individuals"
    )
plot!(sol, vars = collect(N+1:2*N), label="")

u = sol(tsmc, idxs=collect(N+1:2*N))
us = hcat(u.u...)'
plot(tsmc, us-mean_i, label="")

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
