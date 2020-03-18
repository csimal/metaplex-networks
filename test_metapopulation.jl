using LightGraphs
using LightGraphs.SimpleGraphs
using DifferentialEquations
using Plots

include("SI\\metapopulation.jl")

N = 10
P = 1000
volume = 1
frac_infected = 0.01
ninfected = Int(frac_infected*P)
g = complete_graph(N)
g = random_configuration_model(N,fill(3,N))
g = path_graph(N)
V = fill(volume, N)
β = 0.3
μ = 0.1
tmax = 1.0

mp = Metapopulation(g,V,β,μ)
s0 = fill(P, N)
i0 = fill(ninfected, N)
i0 = zeros(Int, N)
s0[1] = P-ninfected
i0[1] = ninfected
Ptot = sum(s0) + sum(i0)

t, s, i = metapopulation_gillespie(mp,s0,i0, nmax = 1000000, tmax=tmax)

plot(t, sum(i, dims=2)/Ptot,
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )

plot(t,i/P,
    line=:steppre,
    label="",
    xlabel = "Time",
    ylabel = "Normalized #infected individuals"
    )

plot(t,s/Ptot, line=:steppre)


t_mc, s_mc, i_mc = metapopulation_montecarlo(mp,s0,i0,tmax=tmax, nmax=150000, nbins=200)

t_mf, s_mf, i_mf, sol = metapopulation_ode(mp, s0, i0, tmax=tmax, saveat=t_mc)

plot(t_mf, i_mf,
    label="",
    xlabel="Time",
    ylabel="#infected individuals",
    linestyle=:dash
    )
plot!(t, i, line=:steppre, label="")

plot(t_mc, i_mc,
    label="",
    xlabel="Time",
    ylabel="#infected individuals"
    )
plot!(t_mf, i_mf, label="", linestyle=:dash)

plot(t_mc, i_mf-i_mc, label="")

plot(t_mc, sum(i_mc, dims=2)/Ptot,
    label="Average",
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    legend=:bottomright
    )
plot!(t_mf, sum(i_mf, dims=2)/Ptot,
    label="Mean Field")


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
