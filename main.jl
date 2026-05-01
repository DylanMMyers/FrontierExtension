include("src/helpers.jl")
include("src/solvers.jl")
include("src/visuals.jl")
include("src/fe.jl")
using Gurobi
using MAT
using SparseArrays
using Plots
gurobi_env = Gurobi.Env()

"""
Parses the flags passed in the command line

Input: cmd line args (if any)

Output: args dict
"""
function parse_args()
    args = Dict(
        "graph" => "graphs/small/Karate.mat",
        "epsi" => 1.0,
        "lflag" => true,
        "prune" => true,
        "constraint" => false,
        "show" => false,
    )

    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--graph"
            args["graph"] = ARGS[i+1]
            i += 2
        elseif ARGS[i] == "--epsi"
            args["epsi"] = parse(Float64, ARGS[i+1])
            i += 2
        elseif ARGS[i] == "--no-prune"
            args["prune"] = false
            args["lflag"] = false;
            i += 1
        elseif ARGS[i] == "--wedge-constraint"
            args["constraint"] = true
            i += 1
        elseif ARGS[i] == "--show-clusters"
            args["show"] = true
            i += 1
        else
            println("Unknown arg: $(ARGS[i])")
            i += 1
        end
    end

    return args
end

"""
Loads graph, assumes .mat format

Input: path to graph

Output: A (adjacency matrix), xy (coordinates of each node for visualization)
"""
function load_graph(path)
    M = matread(path)
    
    if haskey(M, "A")
        A = sparse(M["A"])
    elseif haskey(M, "B")
        A = sparse(M["B"])
    elseif haskey(M, "Problem")
        A = sparse(M["Problem"]["A"])
    else
        error("Unknown matrix key in $path. Keys: $(keys(M))")
    end

    n = size(A, 1)
    m = div(sum(A), 2)
    xy = haskey(M, "xy") ? M["xy"] : spectral_xy(A)

    println("Loaded $path: n=$n, m=$m")
    return A, xy
end

"""
Prunes the set of lambdas

Input: FE cover

Output: pruned FE cover (lol)
"""
function prune_cover(C)
    # sort from left endpoint
    sorted = sort(C, by = entry -> entry[3])

    pruned = []
    lam_cursor = 0.0

    while lam_cursor < 1.0
        # find all entries whose left bound covers the cursor
        candidates = filter(entry -> entry[3] <= lam_cursor + 1e-10, sorted)

        if isempty(candidates)
            # if this branch is reached, there is a gap in the cover
            println("Warning: gap in cover at lam = $lam_cursor")
            break
        end

        # pick interval which intersects with our interval and reaches furthest right
        best = sort(candidates, by = entry -> entry[4], rev=true)[1]
        push!(pruned, best)
        lam_cursor = best[4]
    end

    println("Pruned $(length(C) - length(pruned)) solutions (from $(length(C)))")
    return pruned
end

"""
Interface for running FE, contains flags for certain utils user may want (see readme)
"""
function main()
    args = parse_args()

    A, xy = load_graph(args["graph"])
    n = size(A, 1)
    w = ones(n)

    println("\nRunning FE with epsi=$(args["epsi"]), lflag=$(args["lflag"]), cflag=$(args["constraint"])")
    t = @elapsed C = FE(A, w, args["epsi"], args["lflag"], args["constraint"])
    println("FE runtime: $(round(t, digits=2))s")

    println("\nRaw Cover:")
    for (i, entry) in enumerate(C)
        println("  Solution $i: lam=$(round(entry[2], digits=6)), range=[$(round(entry[3], digits=6)), $(round(entry[4], digits=6))]")
    end

    if args["prune"] && args["lflag"]
        println("\nPruning:")
        C = prune_cover(C)
        println("\nPruned Cover:")
        for (i, entry) in enumerate(C)
            println("  Solution $i: lam=$(round(entry[2], digits=6)), range=[$(round(entry[3], digits=6)), $(round(entry[4], digits=6))]")
        end
    elseif args["prune"] && !args["lflag"]
        println("\nSkipping prune: need --lflag to compute left bounds for pruning")
    end

    if args["show"]
        println("\nRounded Clusterings:")
        for (i, entry) in enumerate(C)
            lam = entry[2]
            Elist = entry[5]
            soln = entry[6]
            obj, clus, avg_obj, avg_time = CL_round_LP_new(A, Elist, soln, 50, lam)
            println("Solution $i: λ=$(round(lam, digits=6)), $(length(unique(clus))) clusters, obj=$obj")
        end

        println("\nCreating visualizations:")
        for (i, entry) in enumerate(C)
            lam = entry[2]
            Elist = entry[5]
            soln = entry[6]
            _, clus, _, _ = CL_round_LP_new(A, Elist, soln, 50, lam)
            f = display_graph_clustering(A, xy, clus)
            graphname = splitext(basename(args["graph"]))[1]
            savefig(f, "figures/$(graphname)_clustering_$(i)_lam_$(round(lam, digits=4)).pdf")
            println("Saved figure for solution $i, lam=$(round(lam, digits=4))")
        end
    end
end

main()