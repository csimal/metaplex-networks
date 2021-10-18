using NetworkEpidemics
using LightGraphs
using LinearAlgebra

function linearized_system(mp::HeterogeneousMetapopulation{SIS})
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    D = mp.D
    #w = Vector{Float64}(undef, nv(mp.h))
    L = laplacian_matrix(mp.h)
    eig = eigen(Array(L))
    βs = β * mp.ks / (mp.N/nv(mp.h))
    f! = function(du,u,p,t)
        mul!(du, L, u, -D[2], 0.0)
        @. du += βs * u * (1 - u) - (γ*u)
    end
    return f!
end

function bamboozle(mp::HeterogeneousMetapopulation{SIS}, i)
    β = mp.dynamics.β
    γ = mp.dynamics.γ
    D = mp.D
    L = laplacian_matrix(mp.h)
    eig = eigen(Array(L))
    βs = fill(β, nv(mp.h))
    βs[i] = -γ + eig.values[i] + 5
    f! = function(du, u, p, t)
        mul!(du, L, u, -D[2], 0.0)
        @. du += p[1] * u * (1 - u) - (p[2]*u)
    end
    return f!, (βs, γ)
end