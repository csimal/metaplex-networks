### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 3dd8b742-5d92-493c-a8b4-1240eda1323a
begin
	using Pkg
	Pkg.activate()
end

# ╔═╡ 6b329c80-f68e-11eb-3b74-37bbf987d96e
begin
	using NetworkEpidemics
	using LightGraphs
	using OrdinaryDiffEq
	using Plots
	using Random
	using Statistics
	using LinearAlgebra
	using GraphRecipes
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
    p = sortperm(degree(h))
    #eig = eigen(Array(normalized_laplacian(h)))
	eig = eigen(Array(laplacian_matrix(h)))
    V = Array{Float64,2}(undef, nv(h), nv(h))
    for i in 1:nv(h)
        V[:,i] .= eig.vectors[p,i]
    end
    # for each node, find the eigenvector whose entry for that node is highest in magnitude
    eigenmode = Vector{Int}(undef, nv(h))
    for i in 1:nv(h)
        eigenmode[i] = findmax(abs.(eig.vectors[i,:]))[2]
    end
end

# ╔═╡ c6d532fc-6b83-4526-a9fb-5cbba1862ff2
begin
	cors = [cor(gdistances(h, i), abs.(eig.vectors[:,eigenmode[i]])) for i in 1:nv(h)]
end

# ╔═╡ 40c03939-c179-469e-9e57-e4c7ca788275
begin
	local i = 1
	scatter(abs.(eig.vectors[:,eigenmode[i]]), gdistances(h, i), label="")
end

# ╔═╡ 369cde5f-fd59-4623-9e17-862cf8a0f0a6
begin
	local i = 1
	local ds = gdistances(h, i)
	scatter(abs.(eig.vectors[:,eigenmode[i]]), label="Eigenmode")
	scatter!(ds/maximum(ds), label="Geodesic distance")
end

# ╔═╡ b861304a-5350-4d1f-87ae-98ce3a7f32dc
scatter(degree(h), label="", title="Degrees")

# ╔═╡ 79ba5921-ed0f-434a-aa1f-54b7918c35e8
begin
	histogram(degree(h), label="", title="Degree distribution")
	scatter!(degree_histogram(h), label="")
end

# ╔═╡ f680da2d-bd4b-4969-933a-22f5c8ab76e4
scatter(eig.values, label="", title="Normalized laplacian spectrum")

# ╔═╡ fd60a8ea-33a7-49f4-9968-d64fbbaafc7e
heatmap(abs.(eig.vectors'))

# ╔═╡ 84a578f2-75dd-46c7-9036-f7e0e7b2154d
begin
	local p =heatmap(abs.(V'),
		xlabel="Node, ordered by degree",
		ylabel="Eigenvector, ordered by eigenvalue",
		title = "$pref",
		#color=cgrad(:balance)
	)
	savefig("$(pref)-degree-eigenmode.png")
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

# ╔═╡ 18ff5b7f-978a-4fa7-a872-a34f8d412854
pert = rand(10:60, M)

# ╔═╡ 4608e36d-75b1-40fb-98de-ab32cb6d12cb
begin
	x0 = Array{Int,2}(undef, M, 2)
	#x0[2:M,:] .= [950 50]
	x0[1:M,1] .= 1000
	x0[1:M,1] .-= pert
	x0[1:M,2] .= pert
	#x0[1,:] .= [900, 100]
end

# ╔═╡ 46776763-aca3-4e42-841a-3d891fd4566a
begin
	tmax = 20.0
	nmax = 1000000
	nsims = 100
	nbins = 200
end

# ╔═╡ b5b79925-0f7a-421b-83c7-1d9fa9368c0e
begin
	β = 1.0
	γ = 1.05
	D = 0.8

	critical_k = (N / M) * γ / β
end

# ╔═╡ 5d5b052e-6c60-4e67-bb40-9db51ce1a9c5
begin
	mp = HeterogeneousMetapopulation(h, N, ks, D, SIS(β,γ))
	mp_ = HeterogeneousMetapopulation(h, N, ks_, D, SIS(β,γ))
end

# ╔═╡ b15b77dd-4f09-4c39-9e8d-7d07b2e6b98c
begin
	ts_mf_mp, u_mf_mp = meanfield(mp, x0, tmax=tmax);
	ts_mf_mp_, u_mf_mp_ = meanfield(mp_, x0, tmax=tmax);
	nothing
end

# ╔═╡ c93e4c5f-5b35-4bb1-83b8-0fb9759c5a67
begin
	plot(ts_mf_mp, sum(u_mf_mp[2], dims=2)/N, 
		label="Baseline", 
		legend=:topleft,
		xlabel="time",
		ylabel="fraction of infected population"
	)
	plot!(ts_mf_mp_, sum(u_mf_mp_[2], dims=2)/N, label="Perturbed")
	#savefig("london-transport.pdf")
end

# ╔═╡ ab8a2f99-028c-4086-a2b3-3810f0c14e69
begin
	plot(ts_mf_mp, (eig.vectors * u_mf_mp[2]')', label="")
	#savefig("london-eig-baseline.pdf")
end

# ╔═╡ 376b69db-5f01-45f2-a064-44d8f9f224d7
begin
	plot(ts_mf_mp_, (eig.vectors * u_mf_mp_[2]')', label="")
	#savefig("london-eig-perturbed.pdf")
end

# ╔═╡ 828a2571-bde7-4e59-9d11-49f89fa04337
begin
	plot(ts_mf_mp_, u_mf_mp_[2], label="")
	#savefig("london-nodes-perturbed.pdf")
end

# ╔═╡ 7b04f046-cd8f-45ad-b079-c33c7f5f166c
begin
	v = eig.vectors[:,eigenmode[pertnode]]
	w = u_mf_mp_[2][end,:]
	scatter(abs.(v), 
		label="Eigenvector",
	)
	scatter!(w/maximum(w), label="Final pattern")
end

# ╔═╡ 2dc743d3-31a9-4a04-a9f5-270a62ace4ea
cor(abs.(v), w)

# ╔═╡ 2dcdde5b-e6f1-4f8a-9436-7753bbb3bebe
begin
	scatter(abs.(v), w,
		xscale = :log10,
		yscale = :log10,
		label="")
end

# ╔═╡ d45f928e-0683-4e96-b0b8-85932236f1fe
heatmap(inv(eig.vectors) * Diagonal(ks_) * eig.vectors)

# ╔═╡ 04060122-3e37-43ae-93b9-793a51eab984
idx = 10

# ╔═╡ 70e91e1b-c9be-428c-b544-cc0584360624
degree(h,idx)

# ╔═╡ accc0563-76bb-4d9b-81b5-c1091367dc9e
begin
	local tmax = 25.0
	lin = linearized_system(mp)
	#lin_ = linearized_system(mp_)
	bamb, pars = bamboozle(mp, idx) 
	prob = ODEProblem(lin, x0[:,2]/N, (0.0,tmax))
	#prob_ = ODEProblem(lin_, x0[:,2]/N, (0.0,tmax))
	prob_b = ODEProblem(bamb, x0[:,2]/N, (0.0,300.0), pars)
	sol = solve(prob, Tsit5())
	#sol_ = solve(prob_, Tsit5())
	sol_b = solve(prob_b, Rodas5())
	nothing
end

# ╔═╡ 18ecc8fb-b404-4c1a-8837-83a7004e86c7
sol_b

# ╔═╡ 7733fa0b-dc5b-4b85-a281-dfc7d363c718
f(t, u...) = (t, sum(u)/M)

# ╔═╡ b8ecb346-2237-439c-a6a0-843edb786f54
begin
	local p = plot(sol, vars=(f,0:M...), 
		label="",
		xlabel="Time",
		ylabel = "Fraction of infected population",
		title="Unperturbed system - $pref"
	)
	savefig("$(pref)_unperturbed.png")
	p
end

# ╔═╡ c444c1c4-8379-4fea-9b50-696eede95309
begin
	local p = plot(sol_b, vars=(f,0:M...), 
		label="",
		xlabel="Time",
		ylabel="Fraction of infected population",
		title="Perturbed system - $pref (k = $(degree(h,idx)))"
	)
	savefig("$(pref)_perturbed_$(idx)_d$(degree(h,idx)).png")
	p
end

# ╔═╡ f53f239e-5d66-43d5-9d56-8a7dd270e6d7
begin
	plot(sol_b, label="")
	plot!(sol_b, vars=(0,idx), label="Perturbed Node", color=:red)
end

# ╔═╡ 20283649-51eb-40b7-9880-800142f9bdf3
begin
	local v = abs.(eig.vectors[:,eigenmode[idx]])
	local u = sol_b.u[end]
	local w = abs.(u .- mean(u))
	local neibs = neighbors(h, idx)
	local p = scatter(u, v, label="", 
		xlabel="Final Pattern", 
		ylabel="Eigenvector", 
		legend=:bottomright,
		title = "Final pattern vs. eigenvector ($pref)",
		#xscale = :log10,
		#yscale = :log10,
	)
	scatter!(p, [u[idx]], [v[idx]], label="Perturbed Node")
	scatter!(p, u[neibs], v[neibs], label="Neighbors of perturbed node")
	#plot!(identity, [0.0,maximum(u)], linestyle=:dash, color=:grey, label="")
	savefig("$(pref)_mode_pattern_$(idx).png")
	p
end

# ╔═╡ a8d92931-fe8e-407b-af3a-86a02261f465
function final_infection(mp, i, tmax)
	f, p = bamboozle(mp, i)
	prob = ODEProblem(f, x0[:,2]/N, (0.0, tmax), p)
	sol = solve(prob, Rodas5(), saveat=[tmax])
	return sum(sol.u[end])/M, sol.t[end], sol.retcode
end

# ╔═╡ 591adfd2-1071-4f1d-8a49-b250d97e153f
begin
	us = Vector{Float64}(undef, nv(h))
	ts = Vector{Float64}(undef, nv(h))
	rc = Vector{Symbol}(undef, nv(h))
	for i in 1:nv(h)
		us[i], ts[i], rc[i] = final_infection(mp, i, 500.0)
	end
end

# ╔═╡ 57eb9ba0-29b0-4dda-a95f-edaaef70bcee
ts

# ╔═╡ 21ed4f58-2860-466b-a588-eed29b2640e8
rc

# ╔═╡ 68867d4d-caee-4a9a-9db3-cfc6a112efbe
begin
	local p = scatter(degree(h), us, 
		label="",
		xlabel = "Degree of perturbed node",
		ylabel = "Final Total Infection",
		title = "$pref"
	)
	savefig("$(pref)_deg_vs_pattern.png")
	p
end

# ╔═╡ 70f4ed11-853a-46de-ad26-2b8165f1111b
begin
	local p = scatter(closeness_centrality(h), us, 
		label="",
		xlabel = "Closeness Centrality",
		ylabel = "Final Total Infection",
		title = "$pref"
	)
	savefig("$(pref)_closeness_vs_pattern.png")
	p
end

# ╔═╡ Cell order:
# ╠═3dd8b742-5d92-493c-a8b4-1240eda1323a
# ╠═6b329c80-f68e-11eb-3b74-37bbf987d96e
# ╠═34b6b36a-0795-453d-941e-65d18ba1c837
# ╠═2aaf973e-7ad5-4326-8e77-dc6f2b43fae7
# ╠═42bf98dc-c1fa-4b7c-86b0-12bcbc3e0158
# ╠═1c210948-e88e-4b9c-9ea8-6d18a609fe31
# ╠═c6d532fc-6b83-4526-a9fb-5cbba1862ff2
# ╠═40c03939-c179-469e-9e57-e4c7ca788275
# ╠═369cde5f-fd59-4623-9e17-862cf8a0f0a6
# ╠═b861304a-5350-4d1f-87ae-98ce3a7f32dc
# ╠═79ba5921-ed0f-434a-aa1f-54b7918c35e8
# ╠═f680da2d-bd4b-4969-933a-22f5c8ab76e4
# ╠═fd60a8ea-33a7-49f4-9968-d64fbbaafc7e
# ╠═84a578f2-75dd-46c7-9036-f7e0e7b2154d
# ╠═f3b8f8d8-11db-4ff1-a05b-873c05bd8670
# ╠═c21d0ea3-5317-4eab-a437-99a172caf799
# ╠═5d5b052e-6c60-4e67-bb40-9db51ce1a9c5
# ╠═18ff5b7f-978a-4fa7-a872-a34f8d412854
# ╠═4608e36d-75b1-40fb-98de-ab32cb6d12cb
# ╠═46776763-aca3-4e42-841a-3d891fd4566a
# ╠═b15b77dd-4f09-4c39-9e8d-7d07b2e6b98c
# ╠═b5b79925-0f7a-421b-83c7-1d9fa9368c0e
# ╠═c93e4c5f-5b35-4bb1-83b8-0fb9759c5a67
# ╠═ab8a2f99-028c-4086-a2b3-3810f0c14e69
# ╠═376b69db-5f01-45f2-a064-44d8f9f224d7
# ╠═828a2571-bde7-4e59-9d11-49f89fa04337
# ╠═7b04f046-cd8f-45ad-b079-c33c7f5f166c
# ╠═2dc743d3-31a9-4a04-a9f5-270a62ace4ea
# ╠═2dcdde5b-e6f1-4f8a-9436-7753bbb3bebe
# ╠═d45f928e-0683-4e96-b0b8-85932236f1fe
# ╠═04060122-3e37-43ae-93b9-793a51eab984
# ╠═70e91e1b-c9be-428c-b544-cc0584360624
# ╠═accc0563-76bb-4d9b-81b5-c1091367dc9e
# ╠═18ecc8fb-b404-4c1a-8837-83a7004e86c7
# ╠═7733fa0b-dc5b-4b85-a281-dfc7d363c718
# ╠═b8ecb346-2237-439c-a6a0-843edb786f54
# ╠═c444c1c4-8379-4fea-9b50-696eede95309
# ╠═f53f239e-5d66-43d5-9d56-8a7dd270e6d7
# ╠═20283649-51eb-40b7-9880-800142f9bdf3
# ╠═a8d92931-fe8e-407b-af3a-86a02261f465
# ╠═591adfd2-1071-4f1d-8a49-b250d97e153f
# ╠═57eb9ba0-29b0-4dda-a95f-edaaef70bcee
# ╠═21ed4f58-2860-466b-a588-eed29b2640e8
# ╠═68867d4d-caee-4a9a-9db3-cfc6a112efbe
# ╠═70f4ed11-853a-46de-ad26-2b8165f1111b
