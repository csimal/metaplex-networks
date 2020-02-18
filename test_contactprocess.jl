using LightGraphs
using LightGraphs.SimpleGraphs
using DifferentialEquations
using Plots
using Random
using GraphPlot

include("randutils.jl")
include("contact_process.jl")
include("metapopulation.jl")

# Test that contact process on a complete graph is equivalent to single node metapopulation
N = 100
g = complete_graph(N)
β = 0.3

frac_infected = 0.05
n_infected =  Int(N*frac_infected)

X = falses(N)
X[1:n_infected] .= true # complete graph. choice of nodes makes no difference

ts, Xs, Es = contact_process_gillespie(g, X, β)

plot(ts, sum(Xs, dims=2)/N,
    xlabel="time",
    ylabel="Fraction of infected individuals",
    label="",
    line=:steppre,
    title="Contact Process"
    )

tmean, Xmean, Mmean = contact_process_montecarlo(g, X, β, tmax = 0.30, nbins=400, nsims = 1000)
#Smean = sqrt.(Mmean)
plot(tmean, Xmean/N,
    label="Average (Contact Process)",
    legend=:bottomright,
    xlabel = "time",
    ylabel = "Fraction of infected individuals",
    )

mp = Metapopulation(SimpleGraph(1), [1], β, 0.0)
s0 = [N-n_infected]
i0 = [n_infected]

tsmp, smp, imp = metapopulation_gillespie(mp,s0,i0, nmax = 100, tmax=0.3)

plot(tsmp, sum(imp, dims=2)/N,
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    title="Metapopulation"
    )

tsmp_mean, smp_mean, imp_mean = metapopulation_montecarlo(mp,s0,i0, nmax=100, tmax=0.3, nbins=400)

plot(tmean, Xmean/N,
    label="Average (Contact Process)",
    legend=:bottomright,
    xlabel = "time",
    ylabel = "Fraction of infected individuals",
    )
plot!(tsmp_mean, imp_mean/N,
    label = "Average (Metapopulation)",
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )

plot(tmean, log10.(abs.((Xmean-imp_mean)/N)),
    xlabel = "time",
    ylabel="log10 Deviation")

# Test that reusing the same seed for the rng yields the same result

N = 100
g = erdos_renyi(N, 0.5)
β = 0.3

frac_infected = 0.05
n_infected =  Int(N*frac_infected)

X = falses(N)
X[1:n_infected] .= true # complete graph. choice of nodes makes no difference

Random.seed!(42) # magic number
ts1, Xs1, Es1 = contact_process_gillespie(g, X, β)

Random.seed!(42) # magic number
ts2, Xs2, Es2 = contact_process_gillespie(g, X, β)

ts1 == ts2
Xs1 == Xs2
Es1 == Es2

# Test the influence of network structure on epidemic


# other tests
N = 1000
g = random_configuration_model(N,fill(3,N))
β = 0.3
X = falses(N)
X[rand_combination(N,50)] .= true

ts, Xs, Es = contact_process_gillespie(g, X, β)

plot(ts, sum(Xs, dims=2)/N, xlabel="time", ylabel="Fraction of infected individuals", label="")

gplot(g)
loc_x, loc_y = spring_layout(g)
h = SimpleDiGraph(g)
edgecols = fill(colorant"lightgray", ne(h));
edgecols[Es] .= colorant"red";
gplot(h,
      edgestrokec=edgecols,
      arrowlengthfrac=0.05)

sol = contact_process_ode(g, X, β, tmax=50.0)

plot(sol, label="")
us = hcat(sol.u...)'
plot!(sol.t, sum(us,dims=2)/N, label="Mean Field")

tmean, Xmean, Mmean = contact_process_montecarlo(g, X, β, tmax = 50.0, nbins=400, nsims = 1000)
Smean = sqrt.(Mmean)
plot(tmean, Xmean/N,
    label="Average",
    legend=:bottomright,
    xlabel = "time",
    ylabel = "Fraction of infected individuals"
    )
plot!(tmean, (Xmean+Smean)/N, linestyle=:dash, color=:lightblue, label="Standard Deviation")
plot!(tmean, (Xmean-Smean)/N, linestyle=:dash, color=:lightblue, label="")
plot!(sol.t, sum(us,dims=2)/N, label="Individual Based Approximation")
