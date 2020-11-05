include("utils\\categorical_tree.jl")
include("utils\\randutils.jl")

a = collect(1:5)
ct = CategoricalTree{Int64}(a)
ct[1] = 0
ct
ct[2] = 0
ct
ct[3] = 0
ct
ct[4] = 0
ct
ct[5] = 0
ct
sum(ct)
length(ct)
rand_categorical(ct)
a
N = 1000
counts = zeros(Int, length(ct))
for i in 1:N
    counts[rand_categorical(ct)] += 1
end

using Plots
bar(1:length(ct), counts/N, label="Empirical probability", legend=:topleft)
scatter!(1:length(ct), a/sum(ct), label="theoretical probability")

using Distributions
χ = sum(((counts-N*a/sum(ct)).^2)./(N*a/sum(ct)))
quantile(Chisq(length(a)-1), 0.95)

using Statistics

n = 10000
m = 5000
step = 100
function compare_times(n,m,step)
    N = 2:step:n
    t_a = zeros(length(N))
    t_ct = zeros(length(N))
    sd_a = zeros(length(N))
    sd_ct = zeros(length(N))
    ta = zeros(m)
    tct = zeros(m)
    for i in 1:length(N)
        a = collect(1:N[i])
        ct = CategoricalTree{Int64}(a)
        for k in 1:m
            ta[k] = @elapsed rand_categorical(a)
            tct[k] = @elapsed rand_categorical(ct)
        end
        t_a[i] = mean(ta)
        t_ct[i] = mean(tct)
        sd_a[i] = std(ta)
        sd_ct[i] = std(tct)
    end
    return N, t_a, t_ct, sd_a, sd_ct
end

N, ta, tct, sda, sdct = compare_times(n,m,step)

using Plots
plot(N, log2.(ta), label="Array based algorithm",
    legend=:topleft,
    xlabel="N",
    ylabel="log2 t")
plot!(N, log2.(tct), label="Tree based algorithm")

a = spzeros(5)
a[2] = 2
a[3] = 3
a[5] = 5
a
rand_categorical(a)
collect(zip(findnz(a)))

N = 10000
counts = zeros(Int, length(a))
for i in 1:N
    counts[rand_categorical(a)] += 1
end

using Plots
bar(1:length(a), counts/N, label="Empirical probability", legend=:topleft)
scatter!(1:length(a), a/sum(a), label="Theoretical probability")

using Distributions
χ = sum(((counts-N*a/sum(a)).^2)./(N*a/sum(a)))
quantile(Chisq(length(a)-1), 0.95)
