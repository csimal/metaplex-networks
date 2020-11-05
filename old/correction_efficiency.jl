using NetworkEpidemics
using LightGraphs
using Statistics
using Plots
using ColorSchemes
using ProgressMeter

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SI}) = Metapopulation(mp.h, mp.D, SI(χ*mp.dynamics.β))

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIS}) = Metapopulation(mp.h, mp.D, SIS(χ*mp.dynamics.β, mp.dynamics.γ))

CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIR}) = Metapopulation(mp.h, mp.D, SIR(χ*mp.dynamics.β, mp.dynamics.δ))

function L2_norm(t, u, v)
    sum = 0.0
    dt = 0.0
    n = length(t)
    f = (u-v).^2
    for i in 1:n-1
        dt = t[i+1] - t[i]
        sum += dt*(f[i]+f[i+1])/2
    end
    return sqrt(sum)
end

norm_uniform(t, u, v) = maximum(abs.(u-v))

norm_mean(t, u, v) = mean(abs.(u-v))


function correction_effficiency(mpx::Metaplex, mp::Metapopulation, x0_mpx, x0_mp; tmax=100.0, nmax = 1000, nsims=1000, nbins=100, norm=norm_uniform)
    #ts_av_mpx, u_av_mpx = average(mpx, x0_mpx, tmax=tmax, nmax=nmax, nsims=nsims, nbins=nbins, progressbar=false)
    #ts_av_mp, u_av_mp = average(mp, x0_mp, tmax=tmax, nmax=nmax, nsims=nsims, nbins=nbins, progressbar=false)
    ts = LinRange(0.0,tmax, nbins)
    ts_mf_mpx, u_mf_mpx = meanfield(mpx, x0_mpx, tmax=tmax, saveat=ts)
    ts_mf_mp, u_mf_mp = meanfield(mp, x0_mp, tmax=tmax, saveat=ts)
    #return Linfty_norm(ts_av_mpx, sum(u_av_mpx[2], dims=2), sum(u_av_mp[2], dims=2))/N, Linfty_norm(ts_mf_mpx, sum(u_mf_mpx[2], dims=2), sum(u_mf_mp[2], dims=2))/N
    return norm(ts_mf_mpx, sum(u_mf_mpx[2], dims=2)/N, sum(u_mf_mp[2], dims=2))/N
end

function test(h, N, k, ps, β, D, x0_mpx, x0_mp; tmax=100.0, nmax=1000, nsims=1000, nbins=100, nsamples=50)
    n = length(ps)
    #norm_av = zeros(n)
    norm_mf = zeros(n)
    si = SI(β)
    p = Progress(n*nsamples, dt=5.0)
    for i in 1:n
        for j in 1:nsamples
            g = watts_strogatz(N, k, ps[i])
            χ = mean(degree(g))
            mp = Metapopulation(h, D, SI(χ*β))
            mpx = Metaplex(g, h, D, SI(β))
            mf = correction_effficiency(mpx, mp, x0_mpx, x0_mp, tmax=tmax, nmax=nmax, nsims=nsims, nbins=nbins)
            #norm_av[i] += av
            norm_mf[i] += mf
            ProgressMeter.next!(p)
        end
    end
    #return norm_av/nsamples, norm_mf/nsamples
    return norm_mf/nsamples
end


N = 1000
M = 10
h = complete_graph(M)
β = 0.1
D = 0.1

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = [div(i-1,div(N,M))+1 for i in 1:N]

x0_mp = Array{Int,2}(undef, M,2)
x0_mp[2:M,:] .= [100 0]
x0_mp[1, :] .= [90, 10]

x0_i = fill(1, N)
x0_i[1:10] .= 2
x0_μ = rand(1:M, N)
x0_mp = zeros(Int, M, 2)
for i in 1:N
    x0_mp[x0_μ[i], x0_i[i]] += 1
end
x0_mpx =  [x0_i, x0_μ]

ps = LinRange(0.0, 1.0, 50)
tmax = 50.0
nmax = 5000
nsims = 200
nbins = 500
nsamples = 5

@time norm_mf = test(h, N, 100, ps, β, D, x0_mpx, x0_mp, nmax=nmax, nsims=nsims, nbins=nbins, tmax=tmax, nsamples=nsamples)

scatter(ps, norm_mf,
        label="Mean Field",
        xlabel="p",
        ylabel="Efficiency of mean degree correction")
#scatter!(ps, norm_mf,
#        label="Mean Field")

using DelimitedFiles
open("efficiency.dat", "w") do io
    writedlm(io, [ps norm_av norm_mf])
end

savefig("efficiency.pdf")
