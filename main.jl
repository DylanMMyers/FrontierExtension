include("src/helpers.jl")
include("src/solvers.jl")
include("src/visuals.jl")
include("src/fe.jl")
using Gurobi
gurobi_env = Gurobi.Env() # one env across all files, main acts as an interface

# --- Small test graph (a triangle + a triangle sharing no edges) ---
# Nodes 1-3 form a clique, nodes 4-6 form a clique, no edges between
n = 6
A = sparse(
    [1,1,2,4,4,5,3],
    [2,3,3,5,6,6,4],
    ones(7), n, n
)
A = A + A'

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