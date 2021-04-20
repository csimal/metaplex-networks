using NetworkEpidemics
using LightGraphs
using Statistics
using Plots
using ColorSchemes
using Random

Random.seed!(2021)

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SI}) = Metapopulation(mp.h, mp.D, SI(χ*mp.dynamics.β))

N = 1000
M = 10

h = complete_graph(M)
g = watts_strogatz(N, 100, 0.0)

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
nsims = 200
nbins = 200


ts_mp, output_mp = gillespie(mp, x0_mp, tmax=tmax, nmax=nmax)
ts_mpx, output_mpx = gillespie(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax)

plot(ts_mp, output_mp[2], label="")
plot(ts_mpx, output_mpx[2], label="")

ts = LinRange(0.0, tmax, nbins)
ps = LinRange(0.0, 0.4, 5)

ts_mf_mp, u_mf_mp = meanfield(mp, x0_mp, tmax=tmax, saveat=ts)
ts_mf_mpc, u_mf_mpc = meanfield(mpc, x0_mp, tmax=tmax, saveat=ts)

u_mf_mpx = zeros(length(ps), length(ts))
for i in 1:length(ps), j in 1:100
    g = watts_strogatz(N, 100, ps[i])
    mpx = Metaplex(g, h, D, SI(β))
    t, u = meanfield(mpx, [x0_i,x0_μ], tmax=tmax, saveat=ts)
    u_mf_mpx[i,:] += sum(u[2], dims=2)
end
u_mf_mpx /= 100


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
    plot!(ts, u_mf_mpx[1,:]/N,
    label = "p=0.0",
    color = colors[5]
    )
    plot!(ts, u_mf_mpx[2,:]/N,
    label = "p=0.1",
    color = colors[6]
    )
    plot!(ts, u_mf_mpx[3,:]/N,
    label = "p=0.2",
    color = colors[7]
    )
    plot!(ts, u_mf_mpx[4,:]/N,
    label = "p=0.3",
    color = colors[8]
    )
    plot!(ts, u_mf_mpx[5,:]/N,
    label = "p=0.4",
    color = colors[9]
    )
    plot!(ts, sum(u_mf_mpc[2], dims=2)/N,
    label="",
    color = orange
    )


plot(ts, (sum(u_mf_mpc[2], dims=2) - u_mf_mpx[1,:])/N,
    xlabel="Time",
    ylabel="Residual",
    label = "p=0.0",
    color = colors[5]
    )
    plot!(ts, (sum(u_mf_mpc[2], dims=2) - u_mf_mpx[2,:])/N,
    label = "p=0.1",
    color = colors[6]
    )
    plot!(ts, (sum(u_mf_mpc[2], dims=2) - u_mf_mpx[3,:])/N,
    label = "p=0.2",
    color = colors[7]
    )
    plot!(ts, (sum(u_mf_mpc[2], dims=2) - u_mf_mpx[4,:])/N,
    label = "p=0.3",
    color = colors[8]
    )
    plot!(ts, (sum(u_mf_mpc[2], dims=2) - u_mf_mpx[5,:])/N,
    label = "p=0.4",
    color = colors[9]
    )

cols = map(x->[x.r, x.g, x.b], colors.colors)
or = [orange.r, orange.g, orange.b]

using DelimitedFiles

writedlm("umpc.dat", sum(u_mf_mpc[2], dims=2)/N)
writedlm("umpx.dat", u_mf_mpx/N)
writedlm("ts.dat", ts)
writedlm("colors.dat", [cols..., or])

using MATLAB;


mat"""
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
""";
