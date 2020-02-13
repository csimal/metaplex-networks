using LightGraphs
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Plots
using GraphPlot
using ColorSchemes

function network_automorphism_classes(g::SimpleGraph{<:Integer}; tol=1e-4)
    r = lobpcg(float(adjacency_matrix((g))), true, 1, tol = 1e-10)
    x = copy(r.X[:,1]) # get the dominant eigenvector
    p = sortperm(x)
    x .= x[p]
    k = 1
    classes = [[p[1]]]
    for i in 2:length(x)
        if abs(x[i]-x[i-1]) < tol
            push!(classes[k], p[i])
        else
            k += 1
            push!(classes, [p[i]])
        end
    end
    return classes, r.X[:,1]
end

A = [0 1 1 0 0 0 0 0;
     1 0 1 0 0 0 0 0;
     1 1 0 1 1 1 0 0;
     0 0 1 0 0 0 1 1;
     0 0 1 0 0 0 1 1;
     0 0 1 0 0 0 1 1;
     0 0 0 1 1 1 0 0;
     0 0 0 1 1 1 0 0
    ]

g = SimpleGraph(A)
gplot(g)

cs = network_automorphism_classes(g)

function pairpermutation_matrix(n,i,j)
    P = spzeros(Int64, n,n)
    for k in 1:n
        if k != i || k != j
            P[k,k] = 1
        end
    end
    P[i,j] = 1
    P[j,i] = 1
    return P
end

function checkclasses(g::SimpleGraph{<:Integer}, cs)
    A = adjacency_matrix(g)
    N = nv(g)
    for c in cs
        n = length(c)
        if n > 1
            for i in 2:n # NB incorrect test
                P = pairpermutation_matrix(N,c[1],c[i])
                if A != P*A*P
                    return false
                end
            end
        end
    end
    return true
end

checkclasses(g, cs)

g = graphfamous("karate")
cs, v = network_automorphism_classes(g)
ps = [(k,v[k]) for k in vcat(classes...)]
classes = filter(c->length(c)>1, cs)
colors = ColorSchemes.Set1_4.colors;
nodecols = fill(colorant"lightgray", nv(g));
for i in 1:length(classes), v in classes[i]
    nodecols[v] = colors[i];
end
gplot(g, nodelabel=1:34, nodefillc = nodecols)

function test(N)
    p = 0.75
    fails = []
    for i in 1:100
        h = erdos_renyi(N,p)
        cs2, v = network_automorphism_classes(h)
        if !checkclasses(h,cs2)
            #println("Oh no!")
            push!(fails, (h,cs2,v))
        end
    end
    return fails
end

N = 20
t = test(N)

using GraphPlot

h = t[1][1]
t[1][2]
t[1][3]

gplot(t[1][1], nodelabel = 1:N)
neighbors(h, 8)
neighbors(h, 15)
