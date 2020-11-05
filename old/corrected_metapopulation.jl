using NetworkEpidemics
using LightGraphs
using Statistics
using Plots
using ColorSchemes

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SI}) = Metapopulation(mp.h, mp.D, SI(χ*mp.dynamics.β))

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIS}) = Metapopulation(mp.h, mp.D, SIS(χ*mp.dynamics.β, mp.dynamics.γ))

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIR}) = Metapopulation(mp.h, mp.D, SIR(χ*mp.dynamics.β, mp.dynamics.δ))


N = 1000
M = 10

h = complete_graph(M)
g = erdos_renyi(N, 0.3)
g = random_regular_graph(N, 50)
g = static_scale_free(N, Int(20*N), 2)
g = star_graph(N)
g = watts_strogatz(N, 100, 0.1)
g = cycle_graph(N)

is_connected(g)

k = degree(g)

χ₁ = mean(closeness_centrality(g))
χ₂ = mean(k.^2)/(mean(k)*N)
χ₃ = mean(k)/N #ne(g)/(n*(n-1)/2)
χ₄ = global_clustering_coefficient(g)


β = 0.1
γ = 0.2
D = 0.1

mp = Metapopulation(h, D, SIR(β,γ))
mp_1 = CorrectedMetapopulation(χ₁, mp)
mp_2 = CorrectedMetapopulation(χ₂, mp)
mp_3 = CorrectedMetapopulation(χ₃, mp)
mp_4 = CorrectedMetapopulation(χ₄, mp)

mpx = Metaplex(g, h, D, SIR(β,γ))


x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]

x0_mp = Array{Int,2}(undef, M,3)
x0_mp[2:M,:] .= [100 0 0]
x0_mp[1, :] .= [90, 10, 0]

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = rand(1:M, N)
x0_mp = zeros(Int, M, 3)
for i in 1:N
    x0_mp[x0_μ[i], x0_i[i]] += 1
end

x0_mp

tmax = 50.0
nmax = 1000000
nsims = 200
nbins = 200

ts, output = gillespie(mp, x0_mp, tmax=tmax, nmax=nmax)
ts_mpx, output_mpx = gillespie(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax)

plot(ts, output[2], label="")
plot(ts_mpx, output_mpx[2], label="")

ts_av, u_av = average(mp, x0_mp, tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)
ts_av_1, u_av_1 = average(mp_1, x0_mp, tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)
ts_av_2, u_av_2 = average(mp_2, x0_mp, tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)
ts_av_3, u_av_3 = average(mp_3, x0_mp, tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)
ts_av_4, u_av_4 = average(mp_4, x0_mp, tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)
ts_av_mpx, u_av_mpx = average(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)

ts_mf, u_mf = meanfield(mp, x0_mp, tmax=tmax, saveat=ts_av)
ts_mf_1, u_mf_1 = meanfield(mp_1, x0_mp, tmax=tmax, saveat=ts_av_1)
ts_mf_2, u_mf_2 = meanfield(mp_2, x0_mp, tmax=tmax, saveat=ts_av_2)
ts_mf_3, u_mf_3 = meanfield(mp_3, x0_mp, tmax=tmax, saveat=ts_av_3)
ts_mf_4, u_mf_4 = meanfield(mp_4, x0_mp, tmax=tmax, saveat=ts_av_4)
ts_mf_mpx, u_mf_mpx = meanfield(mpx, [x0_i,x0_μ], tmax=tmax, saveat=ts_av_mpx)

colors = ColorSchemes.tab10;

plot(ts_av, sum(u_av[2], dims=2)/N,
    label="Metapopulation",
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    linestyle=:dash,
    legend=:bottomright,
    fontfamily="Deja Vu",
    color = colors[1]
    )
    plot!(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N,
    label="Metaplex",
    linestyle=:dash,
    color = colors[2]
    )
    plot!(ts_av_1, sum(u_av_1[2], dims=2)/N,
    label="Closeness Centrality",
    linestyle=:dash,
    color = colors[3]
    )
    plot!(ts_av_2, sum(u_av_2[2], dims=2)/N,
    label="Second Moment",
    linestyle=:dash,
    color = colors[4]
    )
    plot!(ts_av_3, sum(u_av_3[2], dims=2)/N,
    label="Mean Degree",
    linestyle=:dash,
    color = colors[5]
    )
    plot!(ts_av_4, sum(u_av_4[2], dims=2)/N,
    label="Clustering Coefficient",
    linestyle=:dash,
    color = colors[6]
    )
    plot!(ts_mf, sum(u_mf[2], dims=2)/N,
    #label="Meanfield (No correction)",
    label = "",
#    xlabel="Time",
#    ylabel="Fraction of infected individuals",
    #legend=:bottomright,
    color = colors[1]
    )
    plot!(ts_mf_mpx, sum(u_mf_mpx[2], dims=2)/N,
#    label="Meanfield (Metaplex)",
    label = "",
    color = colors[2]
    )
    plot!(ts_mf_1, sum(u_mf_1[2], dims=2)/N,
#    label="Meanfield (Closeness Centrality)",
    label = "",
    color = colors[3]
    )
    plot!(ts_mf_2, sum(u_mf_2[2], dims=2)/N,
#    label="Meanfield (Second Moment)",
    label = "",
    color = colors[4]
    )
    plot!(ts_mf_3, sum(u_mf_3[2], dims=2)/N,
#    label="Meanfield (Mean degree)",
    label = "",
    color = colors[5]
    )
    plot!(ts_mf_4, sum(u_mf_4[2], dims=2)/N,
#    label="Meanfield (Clustering Coefficient)",
    label = "",
    color = colors[6]
    )

cols = map(x->[x.r, x.g, x.b], colors.colors)

using MATLAB


mat"""
figure
xlabel('Time')
ylabel('Fraction of Infected Individuals')
hold on
plot($ts_av, $(sum(u_av[2], dims=2)/N), 'LineStyle', '--', 'Color', $(cols[1]))
plot($ts_av_mpx, $(sum(u_av_mpx[2], dims=2)/N), 'LineStyle', '--', 'Color', $(cols[2]))
plot($ts_av_1, $(sum(u_av_1[2], dims=2)/N), 'LineStyle', '--', 'Color', $(cols[3]))
plot($ts_av_2, $(sum(u_av_2[2], dims=2)/N), 'LineStyle', '--', 'Color', $(cols[4]))
plot($ts_av_3, $(sum(u_av_3[2], dims=2)/N), 'LineStyle', '--', 'Color', $(cols[5]))
plot($ts_av_4, $(sum(u_av_4[2], dims=2)/N), 'LineStyle', '--', 'Color', $(cols[6]))

plot($ts_mf, $(sum(u_mf[2], dims=2)/N), 'Color', $(cols[1]))
plot($ts_mf_mpx, $(sum(u_mf_mpx[2], dims=2)/N), 'Color', $(cols[2]))
plot($ts_mf_1, $(sum(u_mf_1[2], dims=2)/N), 'Color', $(cols[3]))
plot($ts_mf_2, $(sum(u_mf_2[2], dims=2)/N), 'Color', $(cols[4]))
plot($ts_mf_3, $(sum(u_mf_3[2], dims=2)/N), 'Color', $(cols[5]))
plot($ts_mf_4, $(sum(u_mf_4[2], dims=2)/N), 'Color', $(cols[6]))
hold off
legend('Metapopulation', 'Metaplex', 'Closeness Centrality', 'Second Moment', 'Mean Degree', 'Clustering Coefficient')
"""
