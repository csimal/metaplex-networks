using NetworkEpidemics
using DifferentialEquations
using Plots
using MATLAB
using Graphs
using Random


include("networks/empirical_networks.jl")
include("linearized.jl")

begin
	Random.seed!(2021)
	#h = contiguous_usa()
	#h = london_transport()
	h = barabasi_albert(100, 5)
	#h = static_scale_free(100, 500, 2)
	#h = erdos_renyi(500, 20/500)
	pref = "BA"
end

begin
	M = nv(h)
	N = M*1000
end

begin
    sorted_nodes = sortperm(Graphs.degree(h))
    #eig = eigen(Array(normalized_laplacian(h)))
	eig = eigen(Array(laplacian_matrix(h)))
    V = Array{Float64,2}(undef, nv(h), nv(h))
    for i in 1:nv(h)
        V[:,i] .= eig.vectors[sorted_nodes,i]
    end
    # for each node, find the eigenvector whose entry for that node is highest in magnitude
    eigenmode = Vector{Int}(undef, M)
    for i in 1:M
        eigenmode[i] = findmax(abs.(eig.vectors[i,:]))[2]
    end
end

ks = fill(300, M)

β = 1.0
γ = 1.05
D = 10.0

mp = HeterogeneousMetapopulation(h, N, ks, D, SIS(β,γ))

uniform_seed = true
idx = 1
seednode = sorted_nodes[idx]

begin
	if uniform_seed
		seed = rand(10:60, M)
	else
		seed = zeros(Int, M)
		seed[seednode] = 50
	end
end

begin
	x0 = Array{Int,2}(undef, M, 2)
	x0[1:M,1] .= 1000
	x0[1:M,1] .-= seed
	x0[1:M,2] .= seed
end;

begin
	local tmax = 400
	lin = linearized_system(mp)
	bamb, pars = bamboozle(mp, seednode) 
	prob = ODEProblem(lin, x0[:,2]/N, (0.0,tmax))
	prob_b = ODEProblem(bamb, x0[:,2]/N, (0.0,tmax), pars)
	sol = solve(prob, Tsit5(), callback=TerminateSteadyState())
	sol_b = solve(prob_b, Rodas5(), callback=TerminateSteadyState())
end;

begin
	local v = abs.(eig.vectors[:,eigenmode[seednode]])
	local u = sol_b.u[end]
	local w = (u .- minimum(u)) / (maximum(u) - minimum(u))
	scatter(v/v[seednode], label="Eigenmode",
		title = "Final pattern vs. eigenvector ($pref), k=$(Graphs.degree(h,seednode))",
		markercolor=:white,
		msize=6,
		msc=:blue)
	scatter!(w, label="Final Pattern")
	scatter!([seednode], [w[seednode]], label="Dense Node")
end

function final_infection(mp, i, tmax)
	f, p = bamboozle(mp, i)
	prob = ODEProblem(f, x0[:,2]/N, (0.0, tmax), p)
	sol = solve(prob, Rodas5(), callback=TerminateSteadyState())
	return sol
end

begin
	us = Vector{Float64}(undef, nv(h))
	ts = Vector{Float64}(undef, nv(h))
	for i in 1:nv(h)
		fi = final_infection(mp, i, 2000.0)
		us[i] = sum(fi.u[end])/M
		ts[i] = fi.t[end]
	end
end

begin
	local p = scatter(Graphs.degree(h), us, 
		label="",
		xlabel = "Degree of perturbed node",
		ylabel = "Final Total Infection",
		title = "$pref"
	)
	#savefig("$(pref)_deg_vs_pattern.png")
	p
end

begin
    ds = Graphs.degree(h)
    cs = closeness_centrality(h)
    mat"""
    figure
    subplot(2,2,1)
    plot($ds, $us,'.')
    xlabel('Degree of modified node')
    ylabel('Final Infection')
    subplot(2,2,2)
    plot($cs, $us,'.')
    xlabel('Closeness centrality of modified node')
    ylabel('Final Infection')
    subplot(2,2,3)
    plot($ds, $ts,'.')
    xlabel('Degree of modified node')
    ylabel('Time to Final Infection')
    subplot(2,2,4)
    plot($cs, $ts,'.')
    xlabel('Closeness centrality of modified node')
    ylabel('Time to Final Infection')
    
    """
end

function final_infection(mp, nodes, ρ, tmax)
	f, params = densify(mp, nodes, ρ)
	prob = ODEProblem(f, x0[:,2]/N, (0.0, tmax), params)
	sol = solve(prob, Rodas5(), callback=TerminateSteadyState())
	return sol
end

low_first(k) = sorted_nodes[1:k]
high_first(k) = sorted_nodes[(M-k+1):M]
random(k) = randperm(nv(h))[1:k]

d = 5 # how much denser are dense nodes
tmax = 2000

begin
	us_low = Vector{Float64}(undef, nv(h))
	ts_low = Vector{Float64}(undef, nv(h))
	us_high = Vector{Float64}(undef, nv(h))
	ts_high = Vector{Float64}(undef, nv(h))
	for i in 1:nv(h)
		fi_low = final_infection(mp, low_first(i), d, tmax)
		us_low[i] = sum(fi_low.u[end])/M
		ts_low[i] = fi_low.t[end]
		fi_high = final_infection(mp, high_first(i), d, tmax)
		us_high[i] = sum(fi_high.u[end])/M
		ts_high[i] = fi_high.t[end]
	end
end

begin
	local x = (1:M)/M
	local p = scatter(x, us_low,
		label = "Low degrees first",
		xlabel = "Fraction of dense nodes",
		ylabel = "Final infection",
		legend = :bottomright,
		markersize = 2.5
	)
	scatter!(p, x, us_high, label="High degrees first", markersize = 2.5)
end

begin
	local x = (1:M)/M
    mat"""
    figure
    subplot(1,2,1)
    hold on
    plot($x, $us_low, '.')
    plot($x, $us_high, '.')
    hold off
    legend('Low degrees first', 'High degrees first')
    xlabel('Fraction of modified nodes')
    ylabel('Final Infection')
    subplot(1,2,2)
    hold on
    plot($x, $ts_low, '.')
    plot($x, $ts_high, '.')
    hold off
    legend('Low degrees first', 'High degrees first')
    xlabel('Fraction of modified nodes')
    ylabel('Time to Final Infection')
    """
end