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

ts, Xs = contact_process_gillespie(g, X, β)

plot(ts, Xs/N,
    xlabel="time",
    ylabel="Fraction of infected individuals",
    label="",
    line=:steppre,
    title="Contact Process"
    )

tmean, Xmean, Mmean = contact_process_montecarlo(g, X, β, tmax = 0.3, nbins=400, nsims = 1000)
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
ts1, Xs1, Es1 = contact_process_gillespie(g, X, β, record=true)

Random.seed!(42) # magic number
ts2, Xs2, Es2 = contact_process_gillespie(g, X, β, record=true)

ts1 == ts2
Xs1 == Xs2
Es1 == Es2

# Test the influence of network structure on epidemic


# other tests
N = 1000
g = complete_graph(N)
g = random_configuration_model(N,fill(3,N))
β = 0.3
X = falses(N)
X[rand_combination(N,50)] .= true

@time ts, Xs = contact_process_gillespie(g, X, β)
@time ts, Xs, Es = contact_process_gillespie(g, X, β, record=true)

plot(ts, Xs/N, xlabel="time", ylabel="Fraction of infected individuals", label="")

sol = contact_process_ode(g, X, β, tmax=20.0)

plot(sol, label="")
us = hcat(sol.u...)'
plot!(sol.t, sum(us,dims=2)/N, label="Mean Field")

@time tmean, Xmean, Mmean = contact_process_montecarlo(g, X, β, tmax = 20.0, nbins=400, nsims = 1000)
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

using Statistics
using ProgressMeter

function compare_methods_time(n,m,step)
    N = 10:step:n
    m_a = zeros(length(N))
    m_s = zeros(length(N))
    m_t = zeros(length(N))
    t_a = zeros(length(N))
    t_s = zeros(length(N))
    t_t = zeros(length(N))
    qt_a = zeros(length(N),4)
    qt_s = zeros(length(N),4)
    qt_t = zeros(length(N),4)
    ta = zeros(m)
    ts = zeros(m)
    tt = zeros(m)
    β = 0.3
    f = 0.05 # at least 5 percent infected
    p = Progress(length(N)*m, dt=1.0)
    for i in 1:length(N)
        g = random_configuration_model(N[i], fill(4,N[i]))
        X = falses(N[i])
        n_inf = N[i] < 20 ? 1 : Int(floor(f*N[i]))
        X[1:n_inf] .= true
        for k in 1:m
            ta[k] = @elapsed contact_process_gillespie(g, X, β, nmax=N[i], method=:array)
            ts[k] = @elapsed contact_process_gillespie(g, X, β, nmax=N[i], method=:sparse)
            tt[k] = @elapsed contact_process_gillespie(g, X, β, nmax=N[i], method=:tree)
            ProgressMeter.next!(p; showvalues= [(:N,N[i]),(:iter,k)])
        end
        m_a[i] = mean(ta)
        m_s[i] = mean(ts)
        m_t[i] = mean(tt)
        t_a[i] = median(ta)
        t_s[i] = median(ts)
        t_t[i] = median(tt)
        qt_a[i,:] = quantile(ta, [0.05,0.25,0.75,0.95])
        qt_s[i,:] = quantile(ts, [0.05,0.25,0.75,0.95])
        qt_t[i,:] = quantile(tt, [0.05,0.25,0.75,0.95])
    end
    return N, m_a, m_s, m_t, t_a, t_s, t_t, qt_a, qt_s, qt_t
end

n = 3000
m = 1000
step = 500

N, ma, ms, mt, ta, ts, tt, qta, qts, qtt = compare_methods_time(n,m,step)

plot(N, ta, label="Median Array", legend=:topleft, linecolor=:blue)
plot!(N, ts, label="Median Sparse", linecolor=:red)
plot!(N, tt, label="Median Tree", linecolor=:green)
plot!(N, ma, label="Mean Array", linecolor=:blue, linestyle=:dash)
plot!(N, ms, label="Mean Sparse", linecolor=:red, linestyle=:dash)
plot!(N, mt, label="Mean Tree", linecolor=:green, linestyle=:dash)
plot!(N, qta, label="",
    linestyle=:dot, linecolor=:blue)
plot!(N, qts, label="",
    linestyle=:dot, linecolor=:red)
plot!(N, qtt, label="",
    linestyle=:dot, linecolor=:green)


scatter(N.^2, ta, label="", color=:blue)
scatter!(N.^2, ts, label="", color=:red)
scatter(N.*log2.(N), tt,
    label="",
    color=:green,
    xlabel="N log(N)",
    ylabel="t(N)")
