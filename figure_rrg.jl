using NetworkEpidemics
using LightGraphs
using Statistics
using Plots
using ColorSchemes
using LaTeXStrings


CorrectedMetapopulation(χ::Real, mp::Metapopulation{SI}) = Metapopulation(mp.h, mp.D, SI(χ*mp.dynamics.β))

##

function rewire_network(g::SimpleGraph, p)
    h = SimpleGraph(g)
    for e in edges(g)
        u = src(e)
        v = dst(e) 
        if rand() < p
            if rand() < 0.5
                w = u
            else
                w = v
            end
            x = rand(setdiff(vertices(h), [u,v]))
            # If there's already an edge, we just swap them, so the graph is unchanged
            if !has_edge(h,w,x)
                rem_edge!(h,u,v)
                add_edge!(h,w,x)
            end
        end
    end
    return h
end

##
N = 1000
M = 10

h = complete_graph(M)
g = watts_strogatz(N, 100, 0.0)
g = random_regular_graph(N, 50)

k = degree(g)

χ = mean(k)/N # doesn't change so only compute it once

β = 0.1
D = 0.1

mp = Metapopulation(h, D, SI(β))
mpc = CorrectedMetapopulation(χ, mp)

mpx = Metaplex(g, h, D, SI(β))

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = rand(1:M, N)
x0_mp = zeros(Int, M, 2)
for i in 1:N
    x0_mp[x0_μ[i], x0_i[i]] += 1
end

tmax = 20.0
nmax = 10000
nsims = 1000
nbins = 200


ts_mp, output_mp = gillespie(mp, x0_mp, tmax=tmax, nmax=nmax)
ts_mpx, output_mpx = gillespie(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax)

plot(ts_mp, output_mp[2], label="")
plot(ts_mpx, output_mpx[2], label="")

ts = LinRange(0.0, tmax, nbins)
ps = LinRange(0.0, 1.0, 5)

ts_av_mpx, u_av_mpx = average(mpx, [x0_i,x0_μ], tmax=tmax, nbins=nbins,nsims=nsims,nmax=nmax)

ts_mf_mp, u_mf_mp = meanfield(mp, x0_mp, tmax=tmax, saveat=ts)
ts_mf_mpc, u_mf_mpc = meanfield(mpc, x0_mp, tmax=tmax, saveat=ts)

u_mf_mpx = Vector{Array{Float64,2}}(undef, length(ps))
d_stds = Vector(undef, length(ps))
for i in 1:length(ps)
    f = rewire_network(g, ps[i])
    d_stds[i] = std(degree(f))
    local mpx = Metaplex(f, h, D, SI(β))
    t, u = meanfield(mpx, [x0_i,x0_μ], tmax=tmax, saveat=ts)
    u_mf_mpx[i] = u[2]
end


colors = ColorSchemes.GnBu_9;
orange = ColorSchemes.Oranges_3[3];

plot(ts, sum(u_mf_mpc[2], dims=2)/N,
    #label="Meanfield (No correction)",
    label = "Mean degree correction",
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    fontfamily="Deja Vu Sans",
    legend=:bottomright,
    color = orange
    )
    plot!(ts, sum(u_mf_mpx[1], dims=2)/N,
    label = "p=0.0",
    color = colors[5]
    )
    plot!(ts, sum(u_mf_mpx[2], dims=2)/N,
    label = "p=0.25",
    color = colors[6]
    )
    plot!(ts, sum(u_mf_mpx[3], dims=2)/N,
    label = "p=0.5",
    color = colors[7]
    )
    plot!(ts, sum(u_mf_mpx[4], dims=2)/N,
    label = "p=0.75",
    color = colors[8]
    )
    plot!(ts, sum(u_mf_mpx[5], dims=2)/N,
    label = "p=1.0",
    color = colors[9]
    )
    plot!(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N,
    label="Metaplex",
    color=orange,
    linestyle=:dash
    )
    plot!(ts, sum(u_mf_mpc[2], dims=2)/N,
    label="",
    color = orange
    )

# necessary for consistent LaTeX
pyplot()

plot(ps, d_stds.^2,
    label="",
    xlabel="p",
    ylabel=L"\sigma^2")

cols = map(x->[x.r, x.g, x.b], colors.colors)
or = [orange.r, orange.g, orange.b]

##
using MATLAB


eval_string("""
figure
xlabel('Time')
ylabel('Fraction of Infected Individuals')
hold on
plot($ts, $(sum(u_mf_mpc[2], dims=2)/N), 'Color', $(or))
plot($ts, $(sum(u_mf_mpx[1], dims=2)/N), 'Color', $(cols[5]))
plot($ts, $(sum(u_mf_mpx[2], dims=2)/N), 'Color', $(cols[6]))
plot($ts, $(sum(u_mf_mpx[3], dims=2)/N), 'Color', $(cols[7]))
plot($ts, $(sum(u_mf_mpx[4], dims=2)/N), 'Color', $(cols[8]))
plot($ts, $(sum(u_mf_mpx[5], dims=2)/N), 'Color', $(cols[9]))
hold off
legend('Mean degree correction', 'p=0.0', 'p=0.1', 'p=0.2', 'p=0.3', 'p=0.4')
""");
