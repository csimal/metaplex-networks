using NetworkEpidemics
using LightGraphs
using Plots

N = 1000
g = complete_graph(N)
g = erdos_renyi(N, 0.3)
g = random_regular_graph(N,15)
M = 10
h = complete_graph(M)

β = 0.1
D = 0.1

# SI

mp = Metapopulation(h, D, SI(β))
mpx = Metaplex(g, h, D, SI(β))

x0_mp = Array{Int,2}(undef, M,2)
x0_mp[2:M,:] .= [100 0]
x0_mp[1, :] .= [90, 10]

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]

tmax = 100.0
nmax = 12000

@time ts_mp, output_mp = gillespie(mp, x0_mp, nmax=nmax)
@time ts_av_mp, u_av_mp = average(mp, x0_mp, nmax=nmax, tmax=tmax, progressbar=true)
@time ts_mf_mp, u_mf_mp  = meanfield(mp, x0_mp, tmax=tmax, saveat=ts_av_mp)

@time ts_mpx, output_mpx = gillespie(mpx, [x0_i,x0_μ], nmax=nmax)
@time ts_av_mpx, u_av_mpx = average(mpx, [x0_i,x0_μ], nmax=nmax, tmax=tmax, progressbar=true)
@time ts_mf_mpx, u_mf_mpx = meanfield(mpx, [x0_i, x0_μ], tmax=tmax, saveat=ts_av_mpx)

plot(ts_mp, output_mp[2],
    line=:steppre,
    label="",
    title="Metapopulation")
plot(ts_mpx, output_mpx[2],
    line=:steppre,
    label="",
    title="Metaplex")

plot(ts_av_mp, u_av_mp[2],
    label="",
    linestyle=:dash,
    title="Metapopulation")
plot!(ts_mf_mp, u_mf_mp[2], label="")

plot(ts_av_mpx, u_av_mpx[2],
    label="",
    linestyle=:dash,
    title="Metaplex")
plot!(ts_mf_mpx, u_mf_mpx[2], label="")

plot(ts_av_mp, sum(u_av_mp[2], dims=2),
    label="Average (Metapopulation)",
    linestyle=:dash,
    legend=:bottomright,
    title="Metapopulation vs Metaplex (SI) - RRG(k=15)"
    )
plot!(ts_mf_mp, sum(u_mf_mp[2], dims=2),
    label="Meanfield (Metapopulation)")
plot!(ts_av_mpx, sum(u_av_mpx[2], dims=2),
    label="Average (Metaplex)",
    linestyle=:dash,
    )
plot!(ts_mf_mpx, sum(u_mf_mpx[2], dims=2),
    label="Meanfield (Metaplex)")


# SIS
β = 0.1
γ = 1.0
D = 0.1

mp = Metapopulation(h, D, SIS(β,γ))
mpx = Metaplex(g, h, D, SIS(β,γ))
c_p = ContactProcess(g, SIS(β,γ))

x0_mp = Array{Int,2}(undef, M,2)
x0_mp[2:M,:] .= [100 0]
x0_mp[1, :] .= [90, 10]

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]

x0_cp = x0_i .== 2

tmax = 10.0
nmax = 8000

@time ts_mp, output_mp = gillespie(mp, x0_mp, nmax=nmax)
@time ts_av_mp, u_av_mp = average(mp, x0_mp, nmax=nmax, tmax=tmax, progressbar=true)
@time ts_mf_mp, u_mf_mp  = meanfield(mp, x0_mp, tmax=tmax, saveat=ts_av_mp)

@time ts_mpx, output_mpx = gillespie(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax)
@time ts_av_mpx, u_av_mpx = average(mpx, [x0_i,x0_μ], nmax=nmax, tmax=tmax, progressbar=true)
@time ts_mf_mpx, u_mf_mpx = meanfield(mpx, [x0_i, x0_μ], tmax=tmax, saveat=ts_av_mpx)

@time ts_cp, output_cp = gillespie(c_p, x0_cp, nmax=nmax)
@time ts_av_cp, u_av_cp = average(c_p, x0_cp, nmax=nmax, tmax=tmax, progressbar=true)
@time ts_mf_cp, u_mf_cp = meanfield(c_p, x0_cp, tmax=tmax, saveat=ts_av_cp)

plot(ts_mp, output_mp[2],
    line=:steppre,
    label="",
    title="Metapopulation")
plot(ts_mpx, output_mpx[2],
    line=:steppre,
    label="",
    title="Metaplex")
plot(ts_cp, output_cp,
    line=:steppre,
    label="",
    title="Contact Process")

plot(ts_av_mp, u_av_mp[2],
    label="",
    linestyle=:dash,
    title="Metapopulation")
plot!(ts_mf_mp, u_mf_mp[2], label="")

plot(ts_av_mpx, u_av_mpx[2],
    label="",
    linestyle=:dash,
    title="Metaplex")
plot!(ts_mf_mpx, u_mf_mpx[2], label="")


plot(ts_av_mp, sum(u_av_mp[2], dims=2)/N,
    label="Average (Metapopulation)",
    linestyle=:dash,
    legend=:bottomright,
    title="Metapopulation vs Metaplex (SIS) - RRG(k=15)",
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )
plot!(ts_mf_mp, sum(u_mf_mp[2], dims=2)/N,
    label="Meanfield (Metapopulation)")
plot!(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N,
    label="Average (Metaplex)",
    linestyle=:dash,
    )
plot!(ts_mf_mpx, sum(u_mf_mpx[2], dims=2)/N,
    label="Meanfield (Metaplex)")
plot!(ts_mf_cp, sum(u_mf_cp, dims=2)/N,
    label="Meanfield (Contact Process)")


# SIR
β = 0.1
γ = 1.0
D = 0.1

mp = Metapopulation(h, D, SIR(β,γ))
mpx = Metaplex(g, h, D, SIR(β,γ))

x0_mp = Array{Int,2}(undef, M, 3)
x0_mp[2:M,:] .= [100 0 0]
x0_mp[1, :] .= [90, 10, 0]

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]

tmax = 18.0
nmax = 5000

@time ts_mp, output_mp = gillespie(mp, x0_mp, nmax=nmax)
@time ts_av_mp, u_av_mp = average(mp, x0_mp, nmax=nmax, tmax=tmax, progressbar=true)
@time ts_mf_mp, u_mf_mp  = meanfield(mp, x0_mp, tmax=tmax, saveat=ts_av_mp)

@time ts_mpx, output_mpx = gillespie(mpx, [x0_i,x0_μ], tmax=tmax, nmax=nmax)
@time ts_av_mpx, u_av_mpx = average(mpx, [x0_i,x0_μ], nmax=nmax, tmax=tmax, progressbar=true)
@time ts_mf_mpx, u_mf_mpx = meanfield(mpx, [x0_i, x0_μ], tmax=tmax, saveat=ts_av_mpx)

plot(ts_mp, output_mp[2],
    line=:steppre,
    label="",
    title="Metapopulation")
plot(ts_mpx, output_mpx[2],
    line=:steppre,
    label="",
    title="Metaplex")

plot(ts_av_mp, u_av_mp[2],
    label="",
    linestyle=:dash,
    title="Metapopulation")
plot!(ts_mf_mp, u_mf_mp[2], label="")

plot(ts_av_mpx, u_av_mpx[2],
    label="",
    linestyle=:dash,
    title="Metaplex")
plot!(ts_mf_mpx, u_mf_mpx[2], label="")


plot(ts_av_mp, sum(u_av_mp[2], dims=2)/N,
    label="Average (Metapopulation)",
    linestyle=:dash,
    legend=:topright,
    title="Metapopulation vs Metaplex (SIR) - RRG(k=15)",
    xlabel="Time",
    ylabel="Fraction of infected individuals"
    )
plot!(ts_mf_mp, sum(u_mf_mp[2], dims=2)/N,
    label="Meanfield (Metapopulation)")
plot!(ts_av_mpx, sum(u_av_mpx[2], dims=2)/N,
    label="Average (Metaplex)",
    linestyle=:dash,
    )
plot!(ts_mf_mpx, sum(u_mf_mpx[2], dims=2)/N,
    label="Meanfield (Metaplex)")

plot(ts_av_mp, sum(u_av_mp[3], dims=2)/N,
    label="Average (Metapopulation)",
    linestyle=:dash,
    legend=:bottomright,
    title="Metapopulation vs Metaplex (SIR) - RRG(k=15)",
    xlabel="Time",
    ylabel="Fraction of Recovered individuals"
    )
plot!(ts_mf_mp, sum(u_mf_mp[3], dims=2)/N,
    label="Meanfield (Metapopulation)")
plot!(ts_av_mpx, sum(u_av_mpx[3], dims=2)/N,
    label="Average (Metaplex)",
    linestyle=:dash,
    )
plot!(ts_mf_mpx, sum(u_mf_mpx[3], dims=2)/N,
    label="Meanfield (Metaplex)")
