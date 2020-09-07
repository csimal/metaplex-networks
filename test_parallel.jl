using Distributed

nprocs()
addprocs(4)

@distributed for i=1:5
    println("Foo $i")
end
