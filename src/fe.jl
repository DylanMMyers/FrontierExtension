using Gurobi
using JuMP
using SparseArrays
using LinearAlgebra

"""
File involves the actual fe algorithm and ORLP solver (for theta+ and theta-)

***note that main will handle pruning (and optional rounding) of clusters in main as it provides more freedom to user*** 
as in, there will be a pipeline written in main if you want to naively just run fe and prune/round your set of clusterings

still, the rest of the codebase will be compartmentalized to provide more freedoms (and utility in future projects or other impls)
"""



"""
Frontier Extension algorithm

Input:
A (adjacency matrix), w (node weight vector), epsi (approximation factor)

Output:
C (vector of tuples, each tuple contains four values: x_star, lam_0, lam_minus, lam_plus)
***x_star is soln, lam_l and lam_r represent left and right bounds where solution is a (1 + epsi) approx within these bounds
"""
function FE(A, w, epsi, lflag=false)
    # define size of A and our return set
    n = size(A, 1)
    C = []

    # starting lam_0 (any smaller lam_0 results in one large cluster, trivial case)
    lam_0 = 4.0 / n^2

    while lam_0 < 1.0
        # solve for x* (optimal clustering)
        _, x_star = LazyExactLambdaCCGurobi(A, w, lam_0, false, false, true)

        # compute theta+ and theta- bounds
        theta_plus = ORLP(x_star, A, lam_0, epsi, true)
        lam_plus = lam_0 + theta_plus

        # paper doesnt compute lam_minus so by default we will use lam_0 as lam_minus, but if you want it pass in lflag = true
        lam_minus = lam_0
        if lflag
            theta_minus = ORLP(x_star, A, lam_0, epsi, false)
            lam_minus -= theta_minus
        end

        push!(C, (x_star, lam_0, max(lam_minus, 0.0), min(lam_plus, 1.0)))

        lam_0 = (1 + epsi) * lam_plus
    end

    # so our final value we have entered into soln covers up to lam_plus, but our new lam_0 is lam_plus * (1 + epsi)
    # then we are relying on our new value of lam_0 to cover [lam_plus, 1] however since it is likely >= 1.0 we would terminate
    # we will just compute the final range which spills past 1.0 but covers [lam_plus, 1] 
    # we do this unless our last lam_plus is exactly 1.0 (then its redundant)
    if length(C) > 0
        last_lam_plus = C[end][4]
        if last_lam_plus < 1.0
            _, x_star = LazyExactLambdaCCGurobi(A, w, last_lam_plus, false, false, true)
            push!(C, (x_star, max(last_lam_plus, 0.0), last_lam_plus, 1.0))
        end
    end

    println("FE complete: $(length(C)) solns in cover")
    return C
end



"""
Optimal Range Solver

Given a flag (false == "backward" (theta-), true == "forward" (theta+)), will calculate the endpoint for which a
clustering is a valid (1 + epsi) approximation
Solves two auxillary LPs' to calculate this

Input:
x* (representing opt solution for some lam_0), A (adj matrix), lam_0, epsi, direction flag

Output:
theta (representing the lowest/largest value which we could change lam_0 by and still have it be a (1 + epsi) approx)

"""
function ORLP(x_star, A, lam_0, epsi, forward)
    n = size(A, 1)
    s = forward ? 1 : -1
    npairs = nchoose2 = div(n * (n - 1), 2)

    # idx map for pairs 
    vmap = -ones(Int, n, n)
    idx = 0
    for j = 1:n
        for i = 1:j-1
            idx += 1
            vmap[i,j] = idx
            vmap[j,i] = idx
        end
    end

    # map x* to a vector
    x_vec = zeros(npairs)
    for j = 1:n
        for i = 1:j-1
            x_vec[vmap[i,j]] = x_star[i,j]
        end
    end

    # build c (weight vector), c_ij = 1 - lam_0 if edge, lam_0 if non edge
    c_vec = zeros(npairs)
    for j = 1:n
        for i = 1:j-1
            if A[i,j] > 0
                c_vec[vmap[i,j]] = 1 - lam_0
            else
                c_vec[vmap[i,j]] = -lam_0
            end
        end
    end

    # build perturbation vector: how does the objective change when we shift lambda?
    d_vec = fill(-1.0, npairs)

    # build constraint
    constr_indices = Vector{Vector{Int}}()
    constr_coeffs = Vector{Vector{Float64}}()
    constr_rhs = Float64[]

    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                # push all combos of x_ij + x_ik >= x_jk

                push!(constr_indices, [vmap[i,j], vmap[i,k], vmap[j,k]])
                push!(constr_coeffs, [-1.0, 1.0, 1.0])
                push!(constr_rhs, 0.0)

                push!(constr_indices, [vmap[i,k], vmap[i,j], vmap[j,k]])
                push!(constr_coeffs, [-1.0, 1.0, 1.0])
                push!(constr_rhs, 0.0)

                push!(constr_indices, [vmap[j,k], vmap[i,j], vmap[i,k]])
                push!(constr_coeffs, [-1.0, 1.0, 1.0])
                push!(constr_rhs, 0.0)
            end
        end
    end

    # constraints for x_ij, must be <= 1 (and x_ji >= -1)
    for p = 1:npairs
        push!(constr_indices, [p])
        push!(constr_coeffs, [-1.0])
        push!(constr_rhs, -1.0)
    end

    nconstr = length(constr_rhs)

    # compute known scalar values

    cTx = dot(c_vec, x_vec)
    dTx = dot(d_vec, x_vec)
    b_vec = constr_rhs

    # setup ORLP with Gurobi
    
    nvars = nconstr + 1 # number of y vars + theta
    theta_idx = nvars # 1-idx pos of theta

    # maximize theta or minimize -theta
    obj = zeros(nvars)
    obj[theta_idx] = 1.0

    # constraints say that y and theta >= 0, enforce this now instead of giving gurobi more constraints
    lb = zeros(nvars)

    aptr = Ref{Ptr{Cvoid}}()
    vtypes = repeat(GRB_CONTINUOUS, nvars)
    err = GRBnewmodel(gurobi_env, aptr, "ORLP", nvars, obj, lb, C_NULL, vtypes, C_NULL)
    model = aptr[]
    GRBsetintparam(GRBgetenv(model), "OutputFlag", 0)
    GRBsetintattr(model, "ModelSense", -1)

    try
        # first constraint:
        # A^T*y <= c + s*theta*d 
        for p = 1:npairs
            col_indices = Int32[]
            col_coeffs = Float64[]

            for j = 1:nconstr
                for (idx_pos, var_idx) in enumerate(constr_indices[j])
                    if var_idx == p
                        push!(col_indices, Int32(j - 1))
                        push!(col_coeffs, constr_coeffs[j][idx_pos])
                    end
                end
            end

            push!(col_indices, Int32(theta_idx - 1))
            push!(col_coeffs, -s * d_vec[p])

            rhs = c_vec[p]

            err = GRBaddconstr(model, length(col_indices), col_indices, col_coeffs, GRB_LESS_EQUAL, rhs, C_NULL)
        end

        # second constraint
        # massage the definition of the LP down to something Gurobi can process, just algebra
        approx_indices = Int32[]
        approx_coeffs = Float64[]

        for j = 1:nconstr
            push!(approx_indices, Int32(j - 1))
            push!(approx_coeffs, (1 + epsi) * b_vec[j])
        end

        theta_coeff = s * (epsi * nchoose2 - dTx)
        push!(approx_indices, Int32(theta_idx - 1))
        push!(approx_coeffs, theta_coeff)

        approx_rhs = cTx - epsi * lam_0 * nchoose2

        err = GRBaddconstr(model, length(approx_indices), approx_indices, approx_coeffs, GRB_GREATER_EQUAL, approx_rhs, C_NULL)

        # finally, solve for theta and return the found value

        GRBupdatemodel(model)
        GRBoptimize(model)

        robj = Ref{Float64}(0.0)
        GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, robj)
        theta = robj[]

        return theta

    finally
        GRBfreemodel(model)
    end
end
