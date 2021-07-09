
function random_spatial_network(distances::AbstractMatrix, ks, f, ρ)
    N = size(distances, 1)
    k = mean(ks)
    ps = (ks * ks') .* f.(distances) / (ρ * k)
    return edge_probability_graph(ps)
end

function edge_probability_graph(ps::AbstractMatrix)
    N = size(ps, 1)
    edges = [(i,j) for i in 1:N, j in 1:N if rand() < ps[i,j]]
    return SimpleGraphFromIterator(edges)
end