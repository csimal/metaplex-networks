### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3dd8b742-5d92-493c-a8b4-1240eda1323a
begin
	using Pkg
	Pkg.activate()
end

# ╔═╡ 6b329c80-f68e-11eb-3b74-37bbf987d96e
begin
	using NetworkEpidemics
	using Graphs
	using DifferentialEquations
	using Plots
	using Random
	using Statistics
	using LinearAlgebra
	using GraphRecipes
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ 34b6b36a-0795-453d-941e-65d18ba1c837
begin
	#cd("F:\\Dev\\metaplex-networks/old")
	include("networks/empirical_networks.jl")
	include("linearized.jl")
end

# ╔═╡ 2aaf973e-7ad5-4326-8e77-dc6f2b43fae7
begin
	Random.seed!(2021)
	#h = contiguous_usa()
	#h = london_transport()
	h = barabasi_albert(100, 5)
	#h = static_scale_free(100, 500, 2)
	#h = erdos_renyi(500, 20/500)
	pref = "BA"
end

# ╔═╡ 42bf98dc-c1fa-4b7c-86b0-12bcbc3e0158
begin
	M = nv(h)
	N = M*1000
end

# ╔═╡ 1c210948-e88e-4b9c-9ea8-6d18a609fe31
begin
    sorted_nodes = sortperm(Graphs.degree(h))
    #eig = eigen(Array(normalized_laplacian(h)))
	eig = eigen(Array(laplacian_matrix(h)))
    V = Array{Float64,2}(undef, nv(h), nv(h))
    for i in 1:nv(h)
        V[:,i] .= eig.vectors[sorted_nodes,i]
    end
    # for each node, find the eigenvector whose entry for that node is highest in magnitude
    eigenmode = Vector{Int}(undef, nv(h))
    for i in 1:nv(h)
        eigenmode[i] = findmax(abs.(eig.vectors[i,:]))[2]
    end
end

# ╔═╡ b861304a-5350-4d1f-87ae-98ce3a7f32dc
scatter(Graphs.degree(h), label="", title="Degrees")

# ╔═╡ 79ba5921-ed0f-434a-aa1f-54b7918c35e8
begin
	histogram(Graphs.degree(h), label="", title="Degree distribution")
	scatter!(degree_histogram(h), label="")
end

# ╔═╡ fd60a8ea-33a7-49f4-9968-d64fbbaafc7e
heatmap(abs.(eig.vectors'))

# ╔═╡ 84a578f2-75dd-46c7-9036-f7e0e7b2154d
begin
	local p = heatmap(abs.(V'),
		xlabel="Node, ordered by degree",
		ylabel="Eigenvector, ordered by eigenvalue",
		title = "$pref",
		#color=cgrad(:balance)
	)
	#savefig("$(pref)-degree-eigenmode.png")
	p
end

# ╔═╡ f3b8f8d8-11db-4ff1-a05b-873c05bd8670
pertnode = 5

# ╔═╡ c21d0ea3-5317-4eab-a437-99a172caf799
begin
	ks = fill(300, M)
	ks_ = copy(ks)
	ks_[pertnode] = 1000
end

# ╔═╡ 46776763-aca3-4e42-841a-3d891fd4566a
begin
	tmax = 20.0
	nmax = 1000000
	nsims = 100
	nbins = 200
end

# ╔═╡ ba6076e4-f623-49f0-8f1f-ecfac6daac81
md"""
β: $(@bind β Slider(0.0:0.05:2.0; default=1.0, show_value=true))

γ: $(@bind γ Slider(0.0:0.05:2.0; default=1.05, show_value=true))

D: $(@bind D Slider(0.0:0.05:10.0; default=1.0, show_value=true))
"""

# ╔═╡ 5d5b052e-6c60-4e67-bb40-9db51ce1a9c5
begin
	mp = HeterogeneousMetapopulation(h, N, ks, D, SIS(β,γ))
	mp_ = HeterogeneousMetapopulation(h, N, ks_, D, SIS(β,γ))
end

# ╔═╡ 30e0d12e-8d78-4848-be24-e37493cb166a
md"""
Dense node: $(@bind idx Slider(1:nv(h), show_value=true))

Uniform seed? $(@bind uniform_seed CheckBox(default=true))
"""

# ╔═╡ e6b05fd7-5199-4c8e-8236-20a7fe92b211
seednode = sorted_nodes[idx] # NB p is the sorting permutation by degree

# ╔═╡ 18ff5b7f-978a-4fa7-a872-a34f8d412854
begin
	if uniform_seed
		seed = rand(10:60, M)
	else
		seed = zeros(Int, M)
		seed[seednode] = 50
	end
end

# ╔═╡ 4608e36d-75b1-40fb-98de-ab32cb6d12cb
begin
	x0 = Array{Int,2}(undef, M, 2)
	#x0[2:M,:] .= [950 50]
	x0[1:M,1] .= 1000
	x0[1:M,1] .-= seed
	x0[1:M,2] .= seed
	#x0[1,:] .= [900, 100]
end;

# ╔═╡ 7981a920-6766-4f73-a58c-bd3631ac16b1
seed

# ╔═╡ 70e91e1b-c9be-428c-b544-cc0584360624
Graphs.degree(h,seednode) 

# ╔═╡ accc0563-76bb-4d9b-81b5-c1091367dc9e
begin
	local tmax = 400
	lin = linearized_system(mp)
	#lin_ = linearized_system(mp_)
	bamb, pars = bamboozle(mp, seednode) 
	prob = ODEProblem(lin, x0[:,2]/N, (0.0,tmax))
	#prob_ = ODEProblem(lin_, x0[:,2]/N, (0.0,tmax))
	prob_b = ODEProblem(bamb, x0[:,2]/N, (0.0,tmax), pars)
	sol = solve(prob, Tsit5(), callback=TerminateSteadyState())
	#sol_ = solve(prob_, Tsit5())
	sol_b = solve(prob_b, Rodas5(), callback=TerminateSteadyState())
	nothing
end

# ╔═╡ 815a899e-0818-41d8-a676-413d6cb5f0c0
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

# ╔═╡ 7733fa0b-dc5b-4b85-a281-dfc7d363c718
f(t, u...) = (t, sum(u)/M)

# ╔═╡ f53f239e-5d66-43d5-9d56-8a7dd270e6d7
begin
	plot(sol_b, label="")
	plot!(sol_b, vars=(0,seednode), label="Modified Node", color=:red)
end

# ╔═╡ 20283649-51eb-40b7-9880-800142f9bdf3
begin
	local v = abs.(eig.vectors[:,eigenmode[seednode]])
	local u = sol_b.u[end]
	local w = abs.(u .- mean(u))
	local neibs = neighbors(h, seednode)
	local p = scatter(u, v, label="", 
		xlabel="Final Pattern", 
		ylabel="Eigenvector", 
		legend=:bottomright,
		title = "Final pattern vs. eigenvector ($pref)",
		#xscale = :log10,
		#yscale = :log10,
	)
	scatter!(p, [u[seednode]], [v[seednode]], label="Perturbed Node")
	scatter!(p, u[neibs], v[neibs], label="Neighbors of perturbed node")
	#plot!(identity, [0.0,maximum(u)], linestyle=:dash, color=:grey, label="")
	#savefig("$(pref)_mode_pattern_$(idx).png")
	p
end

# ╔═╡ a8d92931-fe8e-407b-af3a-86a02261f465
function final_infection(mp, i, tmax)
	f, p = bamboozle(mp, i)
	prob = ODEProblem(f, x0[:,2]/N, (0.0, tmax), p)
	sol = solve(prob, Rodas5(), callback=TerminateSteadyState())
	return sol
end

# ╔═╡ fe881b04-e7cb-456c-b0fd-5554c1ef76e4
begin
	low_first(k) = sorted_nodes[1:k]
	high_first(k) = sorted_nodes[(nv(h)-k+1):nv(h)]
	random(k) = randperm(nv(h))[1:k]
end

# ╔═╡ 03db73b9-d321-4126-8b78-9c82153f3983
md"""

Number of dense nodes: $(@bind n_dense Slider(1:nv(h), show_value=true)) $(@bind dense_nodes Select([low_first => "lowest degrees", high_first => "highest degrees", random => "random"]))

Density: $(@bind d Slider(1:10, default=5, show_value=true))

"""

# ╔═╡ 59d2e023-6194-4179-a77d-2b107a613bf5
begin
	local tmax = 2000
	nodes = dense_nodes(n_dense)
	dens, ps = densify(mp, nodes, d) 
	prob_d = ODEProblem(dens, x0[:,2]/N, (0.0,tmax), ps)
	sol_d = solve(prob_d, Rodas5(), callback=TerminateSteadyState())
	nothing
end

# ╔═╡ 61016659-b296-4027-b89f-32d7373136dc
plot(sol_d, label="")

# ╔═╡ 8d3cc058-5766-4be5-9be7-30cc131dfa96
function final_infection(mp, nodes, ρ, tmax)
	f, params = densify(mp, nodes, ρ)
	prob = ODEProblem(f, x0[:,2]/N, (0.0, tmax), params)
	sol = solve(prob, Rodas5(), callback=TerminateSteadyState())
	return sol
end

# ╔═╡ 591adfd2-1071-4f1d-8a49-b250d97e153f
begin
	us = Vector{Float64}(undef, nv(h))
	ts = Vector{Float64}(undef, nv(h))
	for i in 1:nv(h)
		fi = final_infection(mp, i, 2000.0)
		us[i] = sum(fi.u[end])/M
		ts[i] = fi.t[end]
	end
end

# ╔═╡ 68867d4d-caee-4a9a-9db3-cfc6a112efbe
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

# ╔═╡ 284e14fd-b358-450d-acc1-baa03c4f468b
findmin(us)

# ╔═╡ 70f4ed11-853a-46de-ad26-2b8165f1111b
begin
	local p = scatter(closeness_centrality(h), us, 
		label="",
		xlabel = "Closeness Centrality",
		ylabel = "Final Total Infection",
		title = "$pref"
	)
	#savefig("$(pref)_closeness_vs_pattern.png")
	p
end

# ╔═╡ 928bb87d-d9f4-49da-b3c3-f674ca1c6d26
begin
	local p = scatter(Graphs.degree(h), ts, 
		label="",
		xlabel = "Degree of perturbed node",
		ylabel = "Time to final infection",
		title = "$pref"
	)
	#savefig("$(pref)_deg_vs_pattern.png")
	p
end

# ╔═╡ 8f3785f6-263c-4c3b-b81c-4c41493f3197
begin
	local p = scatter(closeness_centrality(h), ts, 
		label="",
		xlabel = "Closeness Centrality",
		ylabel = "Time to final infection",
		title = "$pref"
	)
	#savefig("$(pref)_closeness_vs_pattern.png")
	p
end

# ╔═╡ 6df2ca41-cc69-438a-8603-3ab758a4c174
begin
	local tmax = 2000
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

# ╔═╡ f7075b5f-1ad7-460c-a47c-a9f6b7092259
begin
	local x = (1:nv(h))/nv(h)
	local p = scatter(x, us_low,
		label = "Low degrees first",
		xlabel = "Fraction of dense nodes",
		ylabel = "Final infection",
		legend = :bottomright,
		markersize = 2.5
	)
	scatter!(p, x, us_high, label="High degrees first", markersize = 2.5)
end

# ╔═╡ 940e9684-2c19-4b46-a695-8be9afe420a7
begin
	local x = (1:nv(h))/nv(h)
	local p = scatter(x, ts_low,
		label = "Low degrees first",
		xlabel = "Fraction of dense nodes",
		ylabel = "Time to Final infection",
		markersize = 2.5
	)
	scatter!(p, x, ts_high, label="High degrees first", markersize = 2.5)
end

# ╔═╡ Cell order:
# ╠═3dd8b742-5d92-493c-a8b4-1240eda1323a
# ╠═6b329c80-f68e-11eb-3b74-37bbf987d96e
# ╠═34b6b36a-0795-453d-941e-65d18ba1c837
# ╠═2aaf973e-7ad5-4326-8e77-dc6f2b43fae7
# ╠═42bf98dc-c1fa-4b7c-86b0-12bcbc3e0158
# ╠═1c210948-e88e-4b9c-9ea8-6d18a609fe31
# ╠═b861304a-5350-4d1f-87ae-98ce3a7f32dc
# ╠═79ba5921-ed0f-434a-aa1f-54b7918c35e8
# ╠═fd60a8ea-33a7-49f4-9968-d64fbbaafc7e
# ╠═84a578f2-75dd-46c7-9036-f7e0e7b2154d
# ╠═f3b8f8d8-11db-4ff1-a05b-873c05bd8670
# ╠═c21d0ea3-5317-4eab-a437-99a172caf799
# ╠═5d5b052e-6c60-4e67-bb40-9db51ce1a9c5
# ╠═18ff5b7f-978a-4fa7-a872-a34f8d412854
# ╠═4608e36d-75b1-40fb-98de-ab32cb6d12cb
# ╠═46776763-aca3-4e42-841a-3d891fd4566a
# ╟─ba6076e4-f623-49f0-8f1f-ecfac6daac81
# ╠═30e0d12e-8d78-4848-be24-e37493cb166a
# ╠═7981a920-6766-4f73-a58c-bd3631ac16b1
# ╠═e6b05fd7-5199-4c8e-8236-20a7fe92b211
# ╟─70e91e1b-c9be-428c-b544-cc0584360624
# ╠═accc0563-76bb-4d9b-81b5-c1091367dc9e
# ╠═815a899e-0818-41d8-a676-413d6cb5f0c0
# ╠═7733fa0b-dc5b-4b85-a281-dfc7d363c718
# ╠═f53f239e-5d66-43d5-9d56-8a7dd270e6d7
# ╟─20283649-51eb-40b7-9880-800142f9bdf3
# ╠═a8d92931-fe8e-407b-af3a-86a02261f465
# ╠═591adfd2-1071-4f1d-8a49-b250d97e153f
# ╟─68867d4d-caee-4a9a-9db3-cfc6a112efbe
# ╠═284e14fd-b358-450d-acc1-baa03c4f468b
# ╟─70f4ed11-853a-46de-ad26-2b8165f1111b
# ╟─928bb87d-d9f4-49da-b3c3-f674ca1c6d26
# ╟─8f3785f6-263c-4c3b-b81c-4c41493f3197
# ╠═fe881b04-e7cb-456c-b0fd-5554c1ef76e4
# ╠═03db73b9-d321-4126-8b78-9c82153f3983
# ╠═59d2e023-6194-4179-a77d-2b107a613bf5
# ╠═61016659-b296-4027-b89f-32d7373136dc
# ╠═8d3cc058-5766-4be5-9be7-30cc131dfa96
# ╠═6df2ca41-cc69-438a-8603-3ab758a4c174
# ╠═f7075b5f-1ad7-460c-a47c-a9f6b7092259
# ╠═940e9684-2c19-4b46-a695-8be9afe420a7
