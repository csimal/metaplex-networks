using LightGraphs
using LightGraphs.SimpleGraphs
using DifferentialEquations
using Plots
using GraphPlot

include("utils\\randutils.jl")
include("SI\\contact_process.jl")
include("SI\\metapopulation.jl")

# Test that contact process on a complete graph is equivalent to single node metapopulation
N = 100
g = complete_graph(N)
β = 0.3

frac_infected = 0.05
n_infected =  Int(N*frac_infected)

X = falses(N)
X[1:n_infected] .= true # complete graph. choice of nodes makes no difference

t, Xs = contact_process_gillespie(g, X, β)

plot(t, Xs/N,
    xlabel="time",
    ylabel="Fraction of infected individuals",
    label="",
    line=:steppre,
    title="Contact Process"
    )

t_mc, X_mc, sd_mc = contact_process_montecarlo(g, X, β, tmax = 0.3, nbins=400, nsims = 1000)

plot(t_mc, X_mc/N,
    label="Average (Contact Process)",
    legend=:bottomright,
    xlabel = "time",
    ylabel = "Fraction of infected individuals",
    )

mpp = Metapopulation(SimpleGraph(1), [1], β, 0.0)
s0 = [N-n_infected]
i0 = [n_infected]

t_mpp, s_mpp, i_mpp = metapopulation_gillespie(mpp,s0,i0, nmax = 100, tmax=0.3)

plot(t_mpp, sum(i_mpp, dims=2)/N,
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    title="Metapopulation"
    )

t_mpp_mc, s_mpp_mc, i_mpp_mc = metapopulation_montecarlo(mpp,s0,i0, nmax=100, tmax=0.3, nbins=400)

plot(t_mc, X_mc/N,
    label="Average (Contact Process)",
    legend=:bottomright,
    xlabel = "time",
    ylabel = "Fraction of infected individuals",
    )
plot!(t_mpp_mc, i_mpp_mc/N,
    label = "Average (Metapopulation)",
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )

plot(t_mc, log10.(abs.((X_mc-i_mpp_mc)/N)),
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
t1, Xs1, Es1 = contact_process_gillespie(g, X, β, record=true)

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
tmax=0.05
@time ts, Xs = contact_process_gillespie(g, X, β)
@time ts, Xs, Es = contact_process_gillespie(g, X, β, record=true)

plot(ts, Xs/N, xlabel="time", ylabel="Fraction of infected individuals", label="")



@time t_mc, X_mc, sd_mc = contact_process_montecarlo(g, X, β, tmax = tmax, nbins=400, nsims = 1000)

t_mf, X_mf, sol = contact_process_ode(g, X, β, tmax=0.05, saveat=t_mc)

plot(t_mf, X_mf, label="")

plot(t_mf, sum(X_mf, dims=2)/N, label="Mean Field")

plot(t_mc, X_mc/N,
    label="Average",
    legend=:bottomright,
    xlabel = "time",
    ylabel = "Fraction of infected individuals"
    )
plot!(t_mc, (X_mc+sd_mc)/N, linestyle=:dash, color=:lightblue, label="Standard Deviation")
plot!(t_mc, (X_mc-sd_mc)/N, linestyle=:dash, color=:lightblue, label="")
plot!(t_mf, sum(X_mf,dims=2)/N, label="Individual Based Mean Field")

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
