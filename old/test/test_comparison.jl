using LightGraphs
using Plots
using DelimitedFiles

include("utils\\categorical_tree.jl")
include("SI\\contact_process.jl")

A = readdlm("test3reg.dat", ',', Int, '\n')

g = SimpleGraph(A)
N = nv(g)
β = 1.0
X = falses(N)
X[1] = true

ts, Xs = contact_process_gillespie(g, X, β)

plot(ts, Xs)

function test()
    t = zeros(1000)
    for i in 1:1000
        ts, Xs = contact_process_gillespie(g, X, β, tmax=9.0)
        k = 1
        while Xs[k]/N < 0.5
            k += 1
        end
        t[i] = ts[k]
    end
    return t
end

t = test()

using Statistics
mean(t)

function test2()
    tmean, Xmean, Mmean = contact_process_montecarlo(g, X, β, tmax = 10.0, nbins=1000, nsims = 1000)
    k = 1
    while Xmean[k]/N < 0.5
        k += 1
    end
    return tmean[k]
end

test2()

tmean, Xmean, Mmean = contact_process_montecarlo(g, X, β, tmax = 10.0, nbins=1000, nsims = 1000)

plot(tmean, Xmean/1000,
    label="Average (Contact Process)",
    legend=:bottomright,
    xlabel = "time",
    ylabel = "Fraction of infected individuals",
    )
