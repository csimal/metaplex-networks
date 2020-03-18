using LightGraphs
using LightGraphs.SimpleGraphs
using Plots
using ColorSchemes

include("SI\\metaplex.jl")
include("SI\\metapopulation.jl")
include("SI\\contact_process.jl")

N = 100
M = 10

g = complete_graph(N*M)
g = erdos_renyi(N*M,0.5)
h = erdos_renyi(M, 0.5)
h = complete_graph(M)
h = path_graph(M)

β = 0.3
D = 0.1

tmax = 7.0
nmax = 50000


# Metapopulation
mpp = Metapopulation(h, fill(1, M), β, D)
s0 = fill(N, M)
i0 = fill(0, M)
i0[1] = 10
s0[1] = N-10

t_mpp, s_mpp, i_mpp = metapopulation_gillespie(mpp, s0, i0, nmax=nmax, tmax=tmax)

plot(t_mpp, i_mpp,
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="# infected"
    )

plot(t_mpp, sum(i_mpp, dims=2)/(N*M),
    label="",
    line=:steppre,
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    title ="Metapopulation (single simulation)"
    )

t_mpp_mc, s_mpp_mc, i_mpp_mc, sd_s_mpp_mc, sd_i_mpp_mc = metapopulation_montecarlo(mpp,s0,i0,tmax=tmax, nmax=nmax, nbins=400)

t_mpp_mf, s_mpp_mf, i_mpp_mf, mpp_mf = metapopulation_ode(mpp, s0, i0, tmax=tmax, saveat=t_mpp_mc)

plot(t_mpp_mc, i_mpp_mc,
    label="",
    linestyle=:dash,
    xlabel="Time",
    ylabel="#infected individuals"
)
    plot!(t_mpp_mf, i_mpp_mf,
    label=""
)

plot(t_mpp_mc, sum(i_mpp_mc, dims=2)/(N*M),
    label="Average (Metapopulation)",
    legend=:bottomright,
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    linestyle=:dash
    )
    plot!(t_mpp_mf, sum(i_mpp_mf, dims=2)/(N*M),
    label="Mean Field (Metapopulation)",
    )
plot!(t_mpp_mc, sum(i_mpp_mc+sd_i_mpp_mc, dims=2)/(N*M), label="", color=:lightblue)
    plot!(t_mpp_mc, sum(i_mpp_mc-sd_i_mpp_mc, dims=2)/(N*M), label="", color=:lightblue)


plot(t_mpp_mc, sum(i_mpp_mf-i_mpp_mc, dims=2)/(N*M))

# Metaplex
mpx = Metaplex(g,h)
Xi0 = falses(N*M)
Xi0[1:10] .= true
Xμ0 = Xμ0 = [div(i-1,N)+1 for i in 1:N*M]

t_mpx, pcs = metaplex_gillespie(mpx, Xi0, Xμ0, β, [D, D], nmax=nmax, tmax=tmax)

plot(t_mpx, pcs[:,2,:], label="", line=:steppre, title="Metaplex (Gillespie)")

t_mpx_mc, s_mpx_mc, i_mpx_mc, sd_s_mpx_mc, sd_i_mpx_mc = metaplex_montecarlo(mpx, Xi0, Xμ0, β, [D, D], tmax=tmax, nmax=nmax)

t_mpx_mf, s_mpx_mf, i_mpx_mf, sμ_mpx_mf, iμ_mpx_mf, mpx_mf = metaplex_ode(mpx, Xi0, Xμ0, β, [D, D], tmax=tmax, saveat=t_mpx_mc)

plot(t_mpx_mc, i_mpx_mc,
    title="Metaplex (Mean vs MF)",
    label="",
    xlabel="Time",
    ylabel="#infected individuals",
    linestyle=:dash
)
plot!(t_mpx_mf, iμ_mpx_mf, label="")

plot(t_mpx_mc, sum(i_mpx_mc, dims=2)/(N*M),
    label="Average (Metaplex)",
    xlabel="Time",
    ylabel="Fraction of infected individuals")
    plot!(t_mpx_mf, sum(iμ_mpx_mf, dims=2)/(N*M),
    label="Mean Field (Metaplex)")
plot!(t_mpx_mc, sum(i_mpx_mc+sd_i_mpx_mc, dims=2)/(N*M), label="", color=:lightblue)
    plot!(t_mpx_mf, sum(i_mpx_mc-sd_i_mpx_mc, dims=2)/(N*M), label="", color=:lightblue)

t_cp_mc, X_cp_mc, sd_cp_mc = contact_process_montecarlo(g, Xi0, β, tmax=tmax, nmax=nmax)

t_cp_mf, X_cp_mf = contact_process_ode(g, Xi0, β, tmax=tmax, saveat=t_cp_mc)

# Compare the two
colors = ColorSchemes.tab20;
plot(t_mpp_mc, sum(i_mpp_mc, dims=2)/(N*M),
    label="Average (Metapopulation)",
    legend=:bottomright,
    color = colors[1],
    xlabel="Time",
    ylabel="Fraction of infected individuals",
    linestyle=:dash
    )
    plot!(t_mpp_mf, sum(i_mpp_mf, dims=2)/(N*M),
    label="Mean Field (Metapopulation)",
    color=colors[1]
    )
    plot!(t_mpx_mc, sum(i_mpx_mc, dims=2)/(N*M),
    label="Average (Metaplex)",
    linestyle=:dashdot,
    color=colors[3]
    )
    plot!(t_mpx_mf, sum(iμ_mpx_mf, dims=2)/(N*M),
    label="Mean Field (Metaplex)",
    color=colors[3]
    )
plot!(t_mpp_mc, sum(i_mpp_mc+sd_i_mpp_mc, dims=2)/(N*M), label="", color=colors[2])
    plot!(t_mpp_mc, sum(i_mpp_mc-sd_i_mpp_mc, dims=2)/(N*M), label="", color=colors[2])
    plot!(t_mpx_mc, sum(i_mpx_mc+sd_i_mpx_mc, dims=2)/(N*M), label="", color=colors[4])
    plot!(t_mpx_mf, sum(i_mpx_mc-sd_i_mpx_mc, dims=2)/(N*M), label="", color=colors[4])


plot!(t_cp_mc, sum(X_cp_mc, dims=2)/(N*M),
    label="Average (Contact Process)",
    linestyle=:dash
)
plot!(t_cp_mf, sum(X_cp_mf, dims=2)/(N*M),
    label="Mean Field (Contact Process)"
)
