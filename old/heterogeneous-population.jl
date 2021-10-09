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
	using Plots
	using Random
	using Statistics
	using LinearAlgebra
	using GraphRecipes
end

# ╔═╡ 34b6b36a-0795-453d-941e-65d18ba1c837
begin
	cd("F:\\Dev\\metaplex-networks/old")
	include("networks/empirical_networks.jl")
end

# ╔═╡ f51d66a4-2184-42c4-9a75-77f91bda6927
Random.seed!(2021)

# ╔═╡ 2aaf973e-7ad5-4326-8e77-dc6f2b43fae7
begin
	#h = contiguous_usa()
	#h = london_transport()
	h = barabasi_albert(100, 5)
	#h = erdos_renyi(500, 20/500)
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
	heatmap(abs.(V'),
		xlabel="Node, ordered by degree",
		ylabel="Eigenmode, ordered by eigenvalue",
		#color=cgrad(:balance)
	)
	#savefig("london-degree-eigenmode.pdf")
end

# ╔═╡ c21d0ea3-5317-4eab-a437-99a172caf799
begin
	ks = fill(500, M)
	ks_ = copy(ks)
	ks_[1] = 1000
end

# ╔═╡ 18ff5b7f-978a-4fa7-a872-a34f8d412854
pert = rand(10:60, M)

# ╔═╡ 77be188c-8ada-4bf7-97ea-f5f7f5fbb1ca
begin
	using DelimitedFiles
	open("barabasi-albert-adjmat.dat", "w") do io
		writedlm(io, Array(adjacency_matrix(h)))
	end
	open("inner-degrees.dat", "w") do io
		writedlm(io, ks)
	end
	open("infected.dat", "w") do io
		writedlm(io, pert)
	end
end

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
	tmax = 200.0
	nmax = 1000000
	nsims = 100
	nbins = 200
end

# ╔═╡ b5b79925-0f7a-421b-83c7-1d9fa9368c0e
begin
	β = 0.105
	γ = 0.5
	D = 0.5

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
	v = eig.vectors[:,eigenmode[1]]
	w = u_mf_mp_[2][end,:]
	scatter(abs.(v), 
		label="Eigenmode",
	)
	scatter!(w/norm(w), label="Final pattern")
end

# ╔═╡ 2dc743d3-31a9-4a04-a9f5-270a62ace4ea
cor(abs.(v), w)

# ╔═╡ 2dcdde5b-e6f1-4f8a-9436-7753bbb3bebe
begin
	scatter(abs.(v), w,
		xscale=:log10,
		yscale=:log10,
		label="")
end

# ╔═╡ d45f928e-0683-4e96-b0b8-85932236f1fe
heatmap(inv(eig.vectors) * Diagonal(ks_) * eig.vectors)

# ╔═╡ Cell order:
# ╠═3dd8b742-5d92-493c-a8b4-1240eda1323a
# ╠═6b329c80-f68e-11eb-3b74-37bbf987d96e
# ╠═34b6b36a-0795-453d-941e-65d18ba1c837
# ╠═f51d66a4-2184-42c4-9a75-77f91bda6927
# ╠═2aaf973e-7ad5-4326-8e77-dc6f2b43fae7
# ╠═77be188c-8ada-4bf7-97ea-f5f7f5fbb1ca
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
