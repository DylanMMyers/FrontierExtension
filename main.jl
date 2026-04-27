include("src/helpers.jl")
include("src/solvers.jl")
include("src/visuals.jl")
include("src/fe.jl")
using MAT
using SparseArrays
using Gurobi
gurobi_env = Gurobi.Env() # one env across all files, main acts as an interface

function run_tests()
    # graph_dir = "graphs/small"
    # graphs = filter(f -> endswith(f, ".mat"), readdir(graph_dir))

    # for g in graphs
    #     M = matread(joinpath(graph_dir, g))
    #     A = sparse(M["A"])
    #     n = size(A, 1)
    #     m = div(sum(A), 2)
    #     println("$g: n=$n, m=$m")
    # end

    M = matread("graphs/small/Karate.mat")
    A = M["A"]
    A = sparse(A)  # ensure it's sparse

    # Get layout coordinates if available, otherwise compute spectral
    if haskey(M, "xy")
        xy = M["xy"]
    else
        xy = spectral_xy(A)
    end

    n = size(A, 1)
    w = ones(n)
    lam = 0.3

    # --- Test 1: Lazy + ILP (exact integer solution) ---
    c1, X1 = LazyExactLambdaCCGurobi(A, w, lam, true, false, true)
    obj1 = LamCC_obj(A, c1, lam)
    println("Lazy ILP:  clusters = $c1,  obj = $obj1")

    # --- Test 2: Lazy + LP relaxation ---
    _, X2 = LazyExactLambdaCCGurobi(A, w, lam, false, false, true)
    lp_obj2 = LamCC_LP_obj(A, X2, lam)
    println("Lazy LP:   lp_obj = $lp_obj2")
    println("  X2 unique values: $(sort(unique(round.(X2, digits=4))))")

    # --- Test 3: Non-lazy + ILP ---
    opt3, obj3, Elist3, soln3 = LambdaCC_Gurobi(A; weak=false, LP=false, lam=lam)
    println("Full ILP:  obj = $obj3,  optimal = $opt3")

    # --- Test 4: Non-lazy + LP ---
    opt4, obj4, Elist4, soln4 = LambdaCC_Gurobi(A; weak=false, LP=true, lam=lam)
    println("Full LP:   obj = $obj4,  optimal = $opt4")

    # --- Sanity check: LP obj should be <= ILP obj ---
    println("\nLP bound ≤ ILP? $(lp_obj2 <= obj1 + 1e-6)")

    # ORLP tests

    lam_0 = 0.3
    _, X = LazyExactLambdaCCGurobi(A, w, lam_0, false, false, true)

    theta_plus = ORLP(X, A, lam_0, 0.2, true)
    theta_minus = ORLP(X, A, lam_0, 0.2, false)

    println("λ₀ = $lam_0")
    println("θ⁺ = $theta_plus, so x* is (1+ε)-approx up to λ = $(lam_0 + theta_plus)")
    println("θ⁻ = $theta_minus, so x* is (1+ε)-approx down to λ = $(lam_0 - theta_minus)")

    # FE algo test

    println("\n--- FE Algorithm ---")
    C = FE(A, w, 0.2, true)
    for (i, entry) in enumerate(C)
        println("Solution $i: λ=$(entry[2]), range=[$(entry[3]), $(entry[4])]")
    end
end

run_tests()