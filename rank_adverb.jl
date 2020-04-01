# APL/J style rank adverb operator

function rankop(f::Function, n::Integer)
    return function(A::AbstractArray)
        if n == 0
            f.(A)
        else
            m = ndims(A)
            mapslices(f, A, dims = (m-n+1):m)
        end
    end
end

⊙(f,n) = rankop(f, n)

A = rand(3,4,5,6)
rankop(sum, 1)(A)
rankop(sum, 2)(A)
rankop(sum, 3)(A)
# can be used as infix
(sum ⊙ 1)(A)
