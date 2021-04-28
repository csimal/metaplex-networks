using NetworkEpidemics
using LightGraphs
using Statistics
using Plots
using ColorSchemes
using Random

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SI}) = Metapopulation(mp.h, mp.D, SI(χ*mp.dynamics.β))
CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIS}) = Metapopulation(mp.h, mp.D, SIS(χ*mp.dynamics.β, mp.dynamics.γ))
CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIR}) = Metapopulation(mp.h, mp.D, SIR(χ*mp.dynamics.β, mp.dynamics.δ))

N = 1000
M = 10

Random.seed!(2020) # for reproducibility

h = erdos_renyi(M, 0.3)
g = static_scale_free(N, Int(20*N), 2)

k = degree(g)

χ₁ = mean(closeness_centrality(g))
χ₂ = mean(k.^2)/(mean(k)*N)
χ₃ = mean(k)/N #ne(g)/(n*(n-1)/2)
χ₄ = global_clustering_coefficient(g)

β = 0.1
γ = 0.3
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

x0_mp = Array{Int,2}(undef, M,2)
x0_mp[2:M,:] .= [100 0]
x0_mp[1, :] .= [90, 10]

# infected seed randomly dispersed over meta-nodes
x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = rand(1:M, N)
x0_mp = zeros(Int, M, 3)
for i in 1:N
    x0_mp[x0_μ[i], x0_i[i]] += 1
end

tmax = 60.0
nmax = 100000
nsims = 1000
nbins = 200

ts, output = gillespie(mp, x0_mp, tmax=tmax, nmax=nmax)
ts_mpx, output_mpx = gillespie(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax)

plot(ts, output[2], label="")
plot(ts_mpx, output_mpx[2], label="")

Random.seed!(2020)

#ts_av, u_av = average(mp, x0_mp, tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)
ts_av_mpx, u_av_mpx = average(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax, nbins=nbins, nsims=nsims)

ts_mf, u_mf = meanfield(mp, x0_mp, tmax=tmax, saveat=ts_av_mpx)
ts_mf_1, u_mf_1 = meanfield(mp_1, x0_mp, tmax=tmax, saveat=ts_av_mpx)
ts_mf_2, u_mf_2 = meanfield(mp_2, x0_mp, tmax=tmax, saveat=ts_av_mpx)
ts_mf_3, u_mf_3 = meanfield(mp_3, x0_mp, tmax=tmax, saveat=ts_av_mpx)
ts_mf_4, u_mf_4 = meanfield(mp_4, x0_mp, tmax=tmax, saveat=ts_av_mpx)
ts_mf_mpx, u_mf_mpx = meanfield(mpx, [x0_i,x0_μ], tmax=tmax, saveat=ts_av_mpx)

colors = ColorSchemes.tab10;

plot(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N + sum(u_av_mpx[3], dims=2)/N,
    label="Metaplex",
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    linestyle=:dash,
    legend=:bottomright,
    color = colors[2]
    )
    plot!(ts_mf, sum(u_mf[2], dims=2)/N + sum(u_mf[3], dims=2)/N,
    #label="Meanfield (No correction)",
    label = "Metapopulation",
#    xlabel="Time",
#    ylabel="Fraction of infected individuals",
    #legend=:bottomright,
    color = colors[1]
    )
    plot!(ts_mf_mpx, sum(u_mf_mpx[2], dims=2)/N + sum(u_mf_mpx[3], dims=2)/N,
#    label="Meanfield (Metaplex)",
    label = "Metaplex IBMF",
    color = colors[2]
    )
    plot!(ts_mf_1, sum(u_mf_1[2], dims=2)/N + sum(u_mf_1[3], dims=2)/N,
#    label="Meanfield (Closeness Centrality)",
    label = "Closeness",
    color = colors[3]
    )
    plot!(ts_mf_2, sum(u_mf_2[2], dims=2)/N + sum(u_mf_2[3], dims=2)/N,
#    label="Meanfield (Second Moment)",
    label = "Second Moment",
    color = colors[4]
    )
    plot!(ts_mf_3, sum(u_mf_3[2], dims=2)/N + sum(u_mf_3[3], dims=2)/N,
#    label="Meanfield (Mean degree)",
    label = "Mean degree",
    color = colors[5]
    )
    plot!(ts_mf_4, sum(u_mf_4[2], dims=2)/N + sum(u_mf_4[3], dims=2)/N,
#    label="Meanfield (Clustering Coefficient)",
    label = "Clustering",
    color = colors[6]
    )

cols = map(x->[x.r, x.g, x.b], colors.colors)

using MATLAB;


mat"""
figure
title('SI')
xlabel('Time')
ylabel('Fraction of Infected Individuals')
hold on
plot($ts_av_mpx, $(sum(u_av_mpx[2], dims=2)/N + sum(u_av_mpx[3], dims=2)/N), 'LineStyle', '--', 'Color', $(cols[2]))

plot($ts_mf, $(sum(u_mf[2], dims=2)/N + sum(u_mf[3], dims=2)/N), 'Color', $(cols[1]))
plot($ts_mf_mpx, $(sum(u_mf_mpx[2], dims=2)/N + sum(u_mf_mpx[3], dims=2)/N), 'Color', $(cols[2]))
plot($ts_mf_1, $(sum(u_mf_1[2], dims=2)/N + sum(u_mf_1[3], dims=2)/N), 'Color', $(cols[3]))
plot($ts_mf_2, $(sum(u_mf_2[2], dims=2)/N + sum(u_mf_2[3], dims=2)/N), 'Color', $(cols[4]))
plot($ts_mf_3, $(sum(u_mf_3[2], dims=2)/N + sum(u_mf_3[3], dims=2)/N), 'Color', $(cols[5]))
plot($ts_mf_4, $(sum(u_mf_4[2], dims=2)/N + sum(u_mf_4[3], dims=2)/N), 'Color', $(cols[6]))
hold off
ylim([0.0,1.0])
legend({'Metaplex','Metapopulation', 'Metaplex IBMF', 'Closeness Centrality', '\$\\langle k^2\\rangle / \\langle k \\rangle\$', '\$\\langle k \\rangle \$', 'Clustering Coefficient'}, 'Interpreter', 'latex', 'Location', 'southeast')
""";
