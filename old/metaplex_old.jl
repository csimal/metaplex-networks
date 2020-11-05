function metaplex_gillespie(mp::Metaplex{T1,T2}, Xi0::BitVector, Xμ0::Vector{Int}, β, D::Vector; nmax=length(Xi0), tmax=100.0, sampling_method=:tree) where {T1,T2}
    N = nv(mp.g) # number of individuals
    g = mp.g
    M = nv(mp.h) # number of populations
    h = mp.h
    #od = map(k-> k>0 ? 1/float(k) : 0, outdegree(h))
    #od = outdegree(h)
    od = outdegree(h) .> 0
    Xi = copy(Xi0) # epidemic state of each individual
    Xμ = copy(Xμ0) # location of each individual
    popcounts = zeros(Int,2,M) # number of of individuals susceptible and infected for each component
    pops_s = [Set{Int}() for i in 1:M] # Population to individuals without traversing the whole list of individuals
    pops_i = [Set{Int}() for i in 1:M]
    for k in 1:N
        if !Xi[k]
            push!(pops_s[Xμ[k]],k)
            #popcounts[1,Xμ[k]] += 1
        else
            push!(pops_i[Xμ[k]],k)
            #popcounts[2,Xμ[k]] += 1
        end
    end
    popcounts[1,:] = map(length, pops_s)
    popcounts[2,:] = map(length, pops_i)
    a = zeros(N+2*M) # aggregate reactions for each individual becoming infected and individuals migrating from each population
    for k in 1:N
        if !Xi[k]
            l = length(filter(i->Xi[i] && Xμ[i]==Xμ[k], inneighbors(g,k)))
            a[k] = β*l # probability of k recieving infection is β time its number of infected incoming neighbors in the same population
            #nreactions += l
        end
    end
    a[N+1:N+M] .= D[1].*popcounts[1,:].*od
    a[N+M+1:N+2*M] .= D[2].*popcounts[2,:].*od
    if sampling_method == :tree
        a = CategoricalTree(a)
    end
    a0 = sum(a)
    ts = zeros(nmax)
    pcs = zeros(Int,nmax,2,M)
    pcs[1,:,:] .= popcounts
    t = 0.0
    n = 1
    while n < nmax && t < tmax && a0 > 0
        τ = log(1/rand())/a0
        k = rand_categorical(a, a0)
        if k in 1:N # node k becomes infected
            Xi[k] = true
            delete!(pops_s[Xμ[k]], k)
            push!(pops_i[Xμ[k]], k)
            popcounts[1,Xμ[k]] = length(pops_s[Xμ[k]])
            popcounts[2,Xμ[k]] = length(pops_i[Xμ[k]])
            a[N+Xμ[k]] = D[1]*popcounts[1,Xμ[k]]*od[Xμ[k]]
            a[N+M+Xμ[k]] += D[2]*popcounts[2,Xμ[k]]*od[Xμ[k]]
            a[k] = 0.0
            for i in Iterators.filter(j->!Xi[j] && Xμ[j]==Xμ[k], outneighbors(g,k))
                a[i] += β
            end
        elseif k in N+1:N+M # a random susceptible node from population k-N moves to another population
            μ = k-N
            ν = rand(outneighbors(h,μ))
            i = rand(pops_s[μ])
            Xμ[i] = ν
            delete!(pops_s[μ], i)
            push!(pops_s[ν], i)
            popcounts[1,μ] = length(pops_s[μ])
            popcounts[1,ν] = length(pops_s[ν])
            a[N+μ] = D[1]*popcounts[1,μ]*od[μ]
            a[N+ν] = D[1]*popcounts[1,ν]*od[ν]
            l = length(filter(j->Xi[j] && Xμ[j]==ν, inneighbors(g,i)))
            a[i] = β*l
        else # a random infected node from k-N-M moves to another population
            μ = k-N-M
            ν = rand(outneighbors(h,μ))
            i = rand(pops_i[μ])
            Xμ[i] = ν
            delete!(pops_i[μ], i)
            push!(pops_i[ν], i)
            for j in Iterators.filter(l->!Xi[l], outneighbors(g,i))
                if Xμ[j] == μ
                    a[j] -= β
                elseif Xμ[j] == ν
                    a[j] += β
                end
            end
            popcounts[2,μ] = length(pops_i[μ])
            popcounts[2,ν] = length(pops_i[ν])
            a[N+M+μ] = D[2]*popcounts[2,μ]*od[μ]
            a[N+M+ν] = D[2]*popcounts[2,ν]*od[ν]
        end
        a0 = sum(a)
        t += τ
        ts[n+1] = t
        pcs[n+1,:,:] .= popcounts
        n += 1
    end
    return ts[1:n], pcs[1:n,:,:]
end
