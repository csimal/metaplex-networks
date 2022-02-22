using Plots
using ColorSchemes
using Graphs
using MATLAB
using NetworkEpidemics
using Parameters: @with_kw


CorrectedMetapopulation(χ::Real, mp::Metapopulation{SI}) = Metapopulation(mp.h, mp.D, SI(χ*mp.dynamics.β))
CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIS}) = Metapopulation(mp.h, mp.D, SIS(χ*mp.dynamics.β, mp.dynamics.γ))
CorrectedMetapopulation(χ::Real, mp::Metapopulation{SIR}) = Metapopulation(mp.h, mp.D, SIR(χ*mp.dynamics.β, mp.dynamics.δ))

@with_kw mutable struct CorrectedMetapopulationFigure{T}
    seed::Int = 2022
    N::Int
    M::Int
    h::SimpleGraph
    g::SimpleGraph
    dynamics::T
    D::Float64
    mp::Metapopulation{T}
    mp_1::Metapopulation{T}
    mp_2::Metapopulation{T}
    mp_3::Metapopulation{T}
    mp_4::Metapopulation{T}
    mpx::Metaplex{T}
    x0_i::Vector{Int}
    x0_μ::Vector{Int}
    x0_mp::Array{Int,2}
    tmax::Float64 = 30.0
    nmax::Int = 100000
    nsims::Int = 1000
    nbins::Int = 200
    ts::Vector{Float64} = []
    u_av_mpx = []
    u_mf = []
    u_mf_1 = []
    u_mf_2 = []
    u_mf_3 = []
    u_mf_4 = []
    u_mf_mpx = []
    colors = ColorSchemes.tab10
end

function CorrectedMetapopulationFigure(g, h, dyn, D, x0_i, x0_μ, x0_mp, kws...)

    mp = Metapopulation(h, D, dyn)
    mp_1 = CorrectedMetapopulation(χ₁, mp)
    mp_2 = CorrectedMetapopulation(χ₂, mp)
    mp_3 = CorrectedMetapopulation(χ₃, mp)
    mp_4 = CorrectedMetapopulation(χ₄, mp)
    mpx = Metaplex(g, h, D, dyn)
    return CorrectedMetapopulationFigure(
        N = nv(g),
        M = nv(h),
        g = g,
        h = h,
        D = D,
        dynamics = dyn,
        x0_i = x0_i,
        x0_μ = x0_μ,
        x0_mp = x0_mp,
        mp = mp,
        mp_1 = mp_1,
        mp_2 = mp_2,
        mp_3 = mp_3,
        mp_4 = mp_4,
        mpx = mpx,
        kws...
    )
end

function compute_figure_data!(cmf::CorrectedMetapopulationFigure)
    ts_av_mpx, u_av_mpx = average(cmf.mpx, [cmf.x0_i,cmf.x0_μ], tmax=cmf.tmax, nmax=cmf.nmax, nbins=cmf.nbins, nsims=cmf.nsims, progressbar=false)
    _, u_mf = meanfield(cmf.mp, cmf.x0_mp, tmax=cmf.tmax, saveat=ts_av_mpx)
    _, u_mf_1 = meanfield(cmf.mp_1, cmf.x0_mp, tmax=cmf.tmax, saveat=ts_av_mpx)
    _, u_mf_2 = meanfield(cmf.mp_2, cmf.x0_mp, tmax=cmf.tmax, saveat=ts_av_mpx)
    _, u_mf_3 = meanfield(cmf.mp_3, cmf.x0_mp, tmax=cmf.tmax, saveat=ts_av_mpx)
    _, u_mf_4 = meanfield(cmf.mp_4, cmf.x0_mp, tmax=cmf.tmax, saveat=ts_av_mpx)
    _, u_mf_mpx = meanfield(cmf.mpx, [cmf.x0_i,cmf.x0_μ], tmax=cmf.tmax, saveat=ts_av_mpx)
    cmf.ts = ts_av_mpx
    cmf.u_av_mpx = u_av_mpx
    cmf.u_mf = u_mf
    cmf.u_mf_1 = u_mf_1
    cmf.u_mf_2 = u_mf_2
    cmf.u_mf_3 = u_mf_3
    cmf.u_mf_4 = u_mf_4
    cmf.u_mf_mpx = u_mf_mpx
end

function make_plot_julia(f::CorrectedMetapopulationFigure)
    p = plot(f.ts, sum(f.u_av_mpx[2], dims=2)/f.N,
        label = "Metaplex",
        xlabel = "Time",
        ylabel = "Fraction of Infected individuals",
        linestyle = :dash,
        legend = :bottomright,
        color = f.colors[2],
        ylims = (0.0,1.0),
        show = false
    )
    plot!(p, f.ts, sum(f.u_mf[2], dims=2)/f.N,
        label = "Metapopulation",
        color = f.colors[1]
    )
    plot!(p, f.ts, sum(f.u_mf_mpx[2], dims=2)/f.N,
    label = "Metaplex IBMF",
    color = f.colors[2]
    )
    plot!(p, f.ts, sum(f.u_mf_1[2], dims=2)/f.N,
    label = "Closeness",
    color = f.colors[3]
    )
    plot!(p, f.ts, sum(f.u_mf_2[2], dims=2)/f.N,
    label = "Second Moment",
    color = f.colors[4]
    )
    plot!(p, f.ts, sum(f.u_mf_3[2], dims=2)/f.N,
    label = "Mean degree",
    color = f.colors[5]
    )
    plot!(p, f.ts, sum(f.u_mf_4[2], dims=2)/f.N,
    label = "Clustering",
    color = f.colors[6]
    )
    return p
end

function make_plot_julia(f::CorrectedMetapopulationFigure{SIR})
    p = plot(f.ts, sum(f.u_av_mpx[2], dims=2)/f.N + sum(f.u_av_mpx[3], dims=2)/f.N,
        label = "Metaplex",
        xlabel = "Time",
        ylabel = "Fraction of Infected individuals",
        linestyle = :dash,
        legend = :bottomright,
        color = f.colors[2],
        ylims = (0.0,1.0),
        show = false
    )
    plot!(p, f.ts, sum(f.u_mf[2], dims=2)/f.N + sum(f.u_mf[3], dims=2)/f.N,
        label = "Metapopulation",
        color = f.colors[1]
    )
    plot!(p, f.ts, sum(f.u_mf_mpx[2], dims=2)/f.N + sum(f.u_mf_mpx[3], dims=2)/f.N,
    label = "Metaplex IBMF",
    color = f.colors[2]
    )
    plot!(p, f.ts, sum(f.u_mf_1[2], dims=2)/f.N + sum(f.u_mf_1[3], dims=2)/f.N,
    label = "Closeness",
    color = f.colors[3]
    )
    plot!(p, f.ts, sum(f.u_mf_2[2], dims=2)/f.N + sum(f.u_mf_2[3], dims=2)/f.N,
    label = "Second Moment",
    color = f.colors[4]
    )
    plot!(p, f.ts, sum(f.u_mf_3[2], dims=2)/f.N + sum(f.u_mf_3[3], dims=2)/f.N,
    label = "Mean degree",
    color = f.colors[5]
    )
    plot!(p, f.ts, sum(f.u_mf_4[2], dims=2)/f.N + sum(f.u_mf_4[3], dims=2)/f.N,
    label = "Clustering",
    color = f.colors[6]
    )
    return p
end

function make_plot_matlab(f::CorrectedMetapopulationFigure)
    ts = f.ts
    N = f.N
    colors = map(x->[x.r,x.g,x.b], f.colors.colors)
    mat"""
    figure
    xlabel('Time')
    ylabel('Fraction of Infected Individuals')
    hold on
    plot($ts, $(sum(f.u_av_mpx[2], dims=2)/N), 'LineStyle', '--', 'Color', $(colors[2]))
    plot($ts, $(sum(f.u_mf[2], dims=2)/N), 'Color', $(colors[1]))
    plot($ts, $(sum(f.u_mf_mpx[2], dims=2)/N), 'Color', $(colors[2]))
    plot($ts, $(sum(f.u_mf_1[2], dims=2)/N), 'Color', $(colors[3]))
    plot($ts, $(sum(f.u_mf_2[2], dims=2)/N), 'Color', $(colors[4]))
    plot($ts, $(sum(f.u_mf_3[2], dims=2)/N), 'Color', $(colors[5]))
    plot($ts, $(sum(f.u_mf_4[2], dims=2)/N), 'Color', $(colors[6]))
    hold off
    ylim([0.0,1.0])
    legend({'Metaplex','Metapopulation', 'Metaplex IBMF', 'Closeness Centrality', '\$\\langle k^2\\rangle / \\langle k \\rangle\$', '\$\\langle k \\rangle \$', 'Clustering Coefficient'}, 'Interpreter', 'latex', 'Location', 'southeast')
    """
end

function make_plot_matlab(f::CorrectedMetapopulationFigure{SIR})
    ts = f.ts
    N = f.N
    colors = map(x->[x.r,x.g,x.b], f.colors.colors)
    mat"""
    figure
    xlabel('Time')
    ylabel('Fraction of Infected Individuals')
    hold on
    plot($ts, $(sum(f.u_av_mpx[2], dims=2)/N + sum(f.u_av_mpx[3], dims=2)/N), 'LineStyle', '--', 'Color', $(colors[2]))
    plot($ts, $(sum(f.u_mf[2], dims=2)/N + sum(f.u_mf[3], dims=2)/N), 'Color', $(colors[1]))
    plot($ts, $(sum(f.u_mf_mpx[2], dims=2)/N + sum(f.u_av_mpx[3], dims=3)/N), 'Color', $(colors[2]))
    plot($ts, $(sum(f.u_mf_1[2], dims=2)/N + sum(f.u_mf_1[3], dims=2)/N), 'Color', $(colors[3]))
    plot($ts, $(sum(f.u_mf_2[2], dims=2)/N + sum(f.u_mf_2[3], dims=2)/N), 'Color', $(colors[4]))
    plot($ts, $(sum(f.u_mf_3[2], dims=2)/N + sum(f.u_mf_3[3], dims=2)/N), 'Color', $(colors[5]))
    plot($ts, $(sum(f.u_mf_4[2], dims=2)/N + sum(f.u_mf_4[3], dims=2)/N), 'Color', $(colors[6]))
    hold off
    ylim([0.0,1.0])
    legend({'Metaplex','Metapopulation', 'Metaplex IBMF', 'Closeness Centrality', '\$\\langle k^2\\rangle / \\langle k \\rangle\$', '\$\\langle k \\rangle \$', 'Clustering Coefficient'}, 'Interpreter', 'latex', 'Location', 'southeast')
    """
end