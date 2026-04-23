using Gurobi
using JuMP
using SparseArrays
using LinearAlgebra

"""
Involves the actual algorithm which returns a set of Lambdas

File handles fe algo itself, pruning, making calls to other solvers/files

Routine is:
init lam_0 = 1/4n^2
use lp solvers to determine x* for lam_0
call ORLP solvers twice to determine theta+ and theta-
add current cluster to set of opt clusters, tie interval to it (e.g. provide something like [x*, l, r])
set lam_0 -> lam_0 + theta+
repeat until lam_0 >= 1
***note that main will handle pruning (and optional rounding) of clusters in main as it provides more freedom to user*** 
"""

