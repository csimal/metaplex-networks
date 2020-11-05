using LightGraphs
using NetworkEpidemics
using Plots
using ColorSchemes

N = 1000
M = 10
β = 0.05
γ = 0.5

h = erdos_renyi(M, 0.5)
g = static_scale_free(N, 50*N, 2)
g = erdos_renyi(N, 0.2)
mp = Metapopulation(h,0.1, SI(β))
mp2 = ContactProcess(g, SI(β))
mp3 = Metaplex(g, h, 0.1, SI(β))


x0 = Array{Int64,2}(undef, M, 2)
x0[1,:] .= [90, 10]
x0[2:M,:] .= [100 0]

x02 = fill(1, N)
x02[1:10] .= 2
x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]

nmax = 5000
tmax = 15.0

t, output = gillespie(mp, x0, nmax=nmax, tmax=tmax)
t2, output2 = gillespie(mp2, [x02,x0_μ], nmax=nmax, tmax=tmax)

t_av, out_av = average(mp, x0, nsims=1000, nbins=200, tmax=tmax, nmax=nmax)
t_av2, out_av2 = average(mp2, x02, nsims=500, nbins=200, tmax=tmax, nmax=nmax)
t_av3, out_av3 = average(mp3, [x02,x0_μ], nsims=500, nbins=200, tmax=tmax, nmax=nmax)

t_mf, out_mf = meanfield(mp, x0, tmax=tmax, saveat=t_av)
t_mf2, out_mf2 = meanfield(mp2, x02, tmax=tmax, saveat=t_av2)
t_mf3, out_mf3 = meanfield(mp3, [x02,x0_μ], tmax=tmax, saveat=t_av2)

colours = ColorSchemes.tab10;

plot(t_mf, sum(out_mf[2], dims=2)/N,
    label="Metapopulation",
    fontfamily="Deja Vu Sans",
    xlabel="Temps",
    ylabel="Fraction d'individus infectes",
    legend=:bottomright,
    color=colours[1]
    )
    plot!(t_mf2, sum(out_mf2, dims=2)/N,
    label="Contact Process",
    color=colours[2])
    plot!(t_mf3, sum(out_mf3[2], dims=2)/N,
    label="Metaplex",
    color=colours[3]
    )
    plot!(t_av, sum(out_av[2], dims=2)/N,
    label="",
    linestyle=:dash,
    color=colours[1]
    )
    plot!(t_av2, out_av2[2]/N,
    label="",
    linestyle=:dash,
    color=colours[2])
    plot!(t_av3, sum(out_av3[2], dims=2)/N,
    label="",
    linestyle=:dash,
    color=colours[3]
    )


savefig("metaplex-si.pdf")
