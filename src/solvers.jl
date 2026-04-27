using Gurobi
using JuMP
using SparseArrays
using LinearAlgebra

"""
This version of LazyExactLambdaCC does not use the JuMP package,
which makes some things faster since it avoids some overhead
"""
function LazyExactLambdaCCGurobi(A,w,lam,exactflag = true,outputflag = true,verbose = true)
    n = size(A,1)
    W = zeros(n,n)
    for i = 1:n
        for j = 1:n
            W[i,j] = A[i,j]-w[i]*w[j]*lam
        end
    end
    #nvars = div(n*(n-1), 2)
    # build varmap and obj
    vmap = -ones(Int,n,n)
    obj = Float64[]
    nvars = 0
    for j=1:n
        for i=1:j-1
            nvars += 1
            vmap[i,j] = nvars-1
            vmap[j,i] = nvars-1
            push!(obj, W[i,j])
        end
    end
    vtypes = repeat(GRB_BINARY, nvars)
    aptr = Ref{Ptr{Cvoid}}()
    err = GRBnewmodel(gurobi_env, aptr, "LazyLambdaCC", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    m = aptr[]

    if outputflag == false
        err = GRBsetintparam(GRBgetenv(m), "OutputFlag", 0)
    end

    try
    
    cind = Int32[0]
    cval = Float64[0]
    if ~exactflag
        for j=1:n
            for i=1:j-1
                if W[i,j] > 0
                    # -x[i,j] <= 0
                    cind[1] = vmap[i,j]
                    cval[1] = -1
                    error = GRBaddconstr(m, 1, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                else
                    # x[i,j] <= 1
                    cind[1] = vmap[i,j]
                    cval[1] = 1
                    error = GRBaddconstr(m, 1, cind, cval, GRB_LESS_EQUAL, 1.0, C_NULL)
                end
            end
        end
    end

    cind = Int32[0,0,0]
    cval = Float64[0,0,0]
    for i = 1:n
        NeighbsI = findall(x->x>0,(A[:,i]))
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    #@constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                    cind[1] = vmap[min(j,k),max(j,k)]
                    cind[2] = vmap[min(i,k),max(i,k)]
                    cind[3] = vmap[min(i,j),max(i,j)]
                    cval[1] = 1
                    cval[2] = -1
                    cval[3] = -1
                    error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                end
            end
        end
    end

    # Find intial first solution
    if verbose
        println("First round of optimization")
    end
    #JuMP.optimize!(m)
    GRBoptimize(m)
    robj = Ref{Float64}(0.0)
    GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)
    obj = robj[]
    soln = zeros(nvars)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
    D = zeros(n,n)
    for j=1:n
        for i=1:j-1
            D[i,j] = soln[vmap[i,j]+1]
            D[j,i] = soln[vmap[i,j]+1]
        end
    end

    while true
        # x will naturally be upper triangular, but 'find_violations'  wants lower triangular
            #D = Matrix(JuMP.value.(x)')

            # Store violating tuples in a vector
            violations = Vector{Tuple{Int,Int,Int}}()

            find_violations!(D,violations)

            # Iterate through and add constraints
            numvi = size(violations,1)
            if verbose
                print("Adding in $numvi violated constraints.")
            end

            # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
            # that needs to be positive in the inequality:
            # x_ab - x_ac - x_bc <= 0.
            # The other two (b and c) could satisfy either b < c or c <= b.
            # We need to make sure we are only placing constraints in the upper
            # triangular portion of the matrix.
            # We do so by just calling min and max on pairs of nodes
            for v in violations
                #assert(v[1]<v[2])
                cind[1] = vmap[v[1],v[2]]
                cind[2] = vmap[min(v[1],v[3]),max(v[1],v[3])]
                cind[3] = vmap[min(v[2],v[3]),max(v[2],v[3])]
                cval[1] = 1
                cval[2] = -1
                cval[3] = -1
                #@show cind, cval
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                #z = sparsevec(cind.+1, cval, nvars)
                #@show z'*soln
                #@constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
            end
            if numvi == 0
                if verbose
                println(" Optimal solution found.")
                end
                break
            end
            if verbose
                println(" And re-solving the optimization problem.")
            end
            GRBupdatemodel(m)
            err = GRBoptimize(m)
            #@show err
            robj = Ref{Float64}(0.0)
            GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)
            obj = robj[]
            #@show "obj after solve", obj
            GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
            for j=1:n
                for i=1:j-1
                    D[i,j] = soln[vmap[i,j]+1]
                    D[j,i] = soln[vmap[i,j]+1]
                end
            end

        end
        if exactflag
            return extractClustering(D), D
        else
            return ones(n), D
        end
    finally
        #@show "freemodel here!"
        GRBfreemodel(m)
    end
end

"""
LambdaCC_Gurobi
(LambdaCC correlation clustering & LambdaSTC)

Use Gurobi WITHOUT JuMP.jl to solve correlation clustering or
LambdaSTC labeling, either optimally, or just the canonical LP relaxation.

This is quite a bit faster than using the JuMP package.

Input:
    A = adjacency matrix (unsigned, unweighted)
    lam = Lambda parameter in (0,1)
Optional parameters:
    weak = false:
        solve correlation clustering objective and not just stc+

    LP = true:
        just solve the LP relaxation, not the exact ILP

    Timelimit: upper limit on how long to let Gurobi run

    Outputflag: whether or not to show Gurobi solver status during solve

Output:
    obj = output objective
    Elist = ordered list of edges in graph
    Evals = variable value for each edge

From the CIKM paper repo
"""
    function LambdaCC_Gurobi(A; weak, outputflag = false, LP, timelimit = Inf,lam)

    n = size(A,1)
    W = ones(n,n)
    for i = 1:n
        for j = 1:n
            if A[i,j] == 0
                W[i,j] = -1*lam
            else
                W[i,j] = (1-lam)
            end
        end
    end
    # build varmap and obj
    vmap = -ones(Int,n,n)
    obj = Float64[]
    nvars = 0
    Elist = zeros(Int64,round(Int64,n*(n-1)/2),2)
    for j=1:n
        for i=1:j-1
            nvars += 1
            Elist[nvars,1] = i
            Elist[nvars,2] = j
            vmap[i,j] = nvars-1
            vmap[j,i] = nvars-1
            push!(obj, W[i,j])
        end
    end
    aptr = Ref{Ptr{Cvoid}}()
    if LP
        vtypes = repeat(GRB_CONTINUOUS, nvars)
        ub = ones(nvars)
        err = GRBnewmodel(gurobi_env, aptr, "CCLP", nvars, obj, C_NULL, ub, vtypes, C_NULL)

    else
        vtypes = repeat(GRB_BINARY, nvars)
        err = GRBnewmodel(gurobi_env, aptr, "ExactCC", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    end

    m = aptr[]
    GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, timelimit)
    GRBsetintparam(GRBgetenv(m), "OutputFlag",outputflag)

    try
        if !weak
            # all triangle inequality constraints
            cind = Int32[0,0,0]
            cval = Float64[0,0,0]
            for i = 1:n
                for j = i+1:n
                    for k = j+1:n
                        cind[1] = vmap[min(i,j),max(i,j)] # xij
                        cind[2] = vmap[min(i,k),max(i,k)] # xik
                        cind[3] = vmap[min(j,k),max(j,k)] # xjk

                        # xij - xik - xjk <= 0.0
                        cval[1] = 1
                        cval[2] = -1
                        cval[3] = -1
                        error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                        # -xij + xik - xjk <= 0.0
                        cval[1] = -1
                        cval[2] = 1
                        cval[3] = -1
                        error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                        # -xij - xik + xjk <= 0.0
                        cval[1] = -1
                        cval[2] = -1
                        cval[3] = 1
                        error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                    end
                end
            end
        else
            T,W = get_wedges_triangles(A)
            # weak constraints for CC
            # so the only triangle inequality contraints we need
            # are x_ij + x_ik >= x_jk if (i,j,k) is an open wedge centered at i
            cind = Int32[0,0,0]
            cval = Float64[0,0,0]
            for t in W
                i = t[1]    # This is always the center
                j = t[2]
                k = t[3]
                cind[1] = vmap[min(i,j),max(i,j)] # xij
                cind[2] = vmap[min(i,k),max(i,k)] # xik
                cind[3] = vmap[min(j,k),max(j,k)] # xjk
                # x[i,j] + x[i,k] >= x[j,k] so xjk - xij - xik <= 0
                cval[1] = -1
                cval[2] = -1
                cval[3] = 1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
            end

        end

    GRBoptimize(m)
    stat = Ref{Int32}(0)
    GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
    #Status codes: https://www.gurobi.com/documentation/9.5/refman/optimization_status_codes.html#sec:StatusCodes
    status = stat[]
    # println("Status = $status")
    if status == 2
        optimal = true
    else
        optimal = false
    end
    robj = Ref{Float64}(0.0)
    GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)
    #     println("Value of robj is: ",robj)
    #     println("Value of nvars is",nvars)

    obj = robj[] + (lam)*(nvars - sum(A)/2)   # adjust by number of non-edges to get the actual LambdaCC objective
    println("Value of obj is ",obj)
    soln = zeros(nvars)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)

    # Elist is a linearization of the node pairs
    # soln[k] is the x variable value for node-par Elist[k]
    # Sanity check if you ever need it for computing obj
    #     obj2 = 0
    #     for k = 1:size(Elist,1)
    #         i = Elist[k,1]
    #         j = Elist[k,2]
    #         if A[i,j] == 0
    #             obj2 += (1-soln[k])*(lam)
    #         else
    #             obj2 += (soln[k])*(1-lam)
    #         end
    #     end
    #     println("Value of object2 is ",obj2)
    return optimal, obj, Elist, soln
    finally
        #@show "freemodel here!"
        GRBfreemodel(m)
    end
end

"""

HELPERS BELOW

"""

"""
# find_violations
# Given a candidate distance matrix D, iterate through all 3-tuples of nodes
# and return the tuples where triangle inequality constraints have been violated.
#
# Output is stored in vector 'violations'
#
# Heuristic speed ups:
#       - Only grab Dij, Dik, Djk once per tuple
#       - if a - b - c > 0, then a>b and a>c. By checking these first
#           we rule out a lot of situations where constraints do not need to be
#           be added, so in all test cases this gave a speed up
#
# Note that we want D to be lower triangular here, though if it is symmetric
# this will work fine as well. We just need to make sure the distance information
# in D is not just stored in the upper triangular portion
"""
function find_violations!(D::Matrix{Float64}, violations::Vector{Tuple{Int,Int,Int}})
  n = size(D,1)

  # We only need this satisfied to within a given tolerance, since the
  # optimization software will only solve it to within a certain tolerance
  # anyways. This can be tweaked if necessary.
  epsi = 1e-8
  @inbounds for i = 1:n-2
       for j = i+1:n-1
          a = D[j,i]
           for k = j+1:n
              b = D[k,i]
              c = D[k,j]
        if a - b > epsi && a - c > epsi && a-b-c > epsi
            push!(violations, (i,j,k))
                # @constraint(m, x[i,j] - x[i,k] - x[j,k] <= 0)
        end
        if b - a > epsi && b - c > epsi && b-a-c > epsi
            push!(violations, (i,k,j))
            # @constraint(m, x[i,k] - x[i,j] - x[j,k] <= 0)
        end

        if c - a > epsi && c-b>epsi && c-a-b > epsi
            push!(violations, (j,k,i))
            # @constraint(m, x[j,k] - x[i,k] - x[i,j] <= 0)
        end
      end
    end
  end
end

# extractClustering
# Given a 0,1 indicator matrix for node clustering, extract the n x 1
# cluster indicator vector
function extractClustering(x)
    # Assuming the triangle inequality results work, we just have to go through
    # each row of x, and if an element hasn't been placed in a cluster yet, it
    # it starts its own new one and adds everything that is with it.
    n = size(x,1)
    NotClustered = fill(true,n)
    c = zeros(n)
    clusnum = 1
    for i = 1:n
        if NotClustered[i]
            for j = i:n
                if x[i,j] < .01 # must be zero, they are clustered together
                    c[j] = clusnum
                    NotClustered[j] = false;
                end
            end
            clusnum +=1
        end
    end
    return round.(Int64,c)
end