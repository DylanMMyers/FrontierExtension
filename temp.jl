using MAT, SparseArrays

for dir in ["graphs/small", "graphs/medium"]
    for f in filter(x -> endswith(x, ".mat"), readdir(dir))
        M = matread(joinpath(dir, f))
        for k in ["A", "B", "Problem"]
            if haskey(M, k)
                A = sparse(k == "Problem" ? M[k]["A"] : M[k])
                n = size(A, 1)
                println("$dir/$f: n=$n, m=$(div(sum(A),2)), avg_d=$(round(sum(A)/n, digits=1))")
                break
            end
        end
    end
end