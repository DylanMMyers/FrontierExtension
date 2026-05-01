using Random
using SparseArrays
using LinearAlgebra

"""
LambdaCC objective check: given a graph A and a clustering c and parameter lam, check the
LambdacC objetive function (the weight of mistakes).

From the CIKM paper repo
"""
function LamCC_obj(A,c,lam)
    n = size(A,1)
    mistakes = 0
    for i = 1:n
        for j = i+1:n
            if A[i,j] == 0 && c[i] == c[j]
                mistakes += lam
            end
            if A[i,j] == 1 && c[i] != c[j]
                mistakes += 1-lam
            end
        end
    end
    return mistakes
end

"""
LambdaCC LP relaxation
"""
function LamCC_LP_obj(A,X,lam)
    n = size(A,1)
    mistakes = 0
    for i = 1:n
        for j = i+1:n
            if A[i,j] == 0 
                mistakes += lam*(1-X[i,j])
            end
            if A[i,j] == 1 
                mistakes += (1-lam)*X[i,j]
            end
        end
    end
    return mistakes
end

"""
The pivot algorithm.

Copied from CIKM paper repo, but likely present in earlier ICML 2022 paper repo or somewhere else.
"""

function permutation_pivot_faster(A)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    Neighbs = ConstructAdj(A,n)
    clusnum = 1

    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return c, round(Int64,Sbin)
end


# From the adjacency matrix, build an adjacency list for the graph
function ConstructAdj(C::SparseMatrixCSC,n::Int64)
    rp = C.rowval
    ci = C.colptr
    Neighbs = Vector{Vector{Int64}}()
    for i = 1:n
        # chop up the rp vector and put it in Neighbs
        push!(Neighbs,rp[ci[i]:ci[i+1]-1])
    end
    return Neighbs

end



"""
Fast version of the pivot algorithm: a random permutation decides the ordering of pivots.
"""
function permutation_pivot(A)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    Neighbs = ConstructAdj(A,n)
    clusnum = 1
    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
            end
        end
        clusnum += 1
    end
    return c
end

"""
Convert an LP solution to an adjacency list (used for rounding fractional lp solutions)
"""
function lp_to_adj(x_star, threshold=0.5)
    n = size(x_star, 1)
    rows = Int[]
    cols = Int[]
    for j = 1:n
        for i = 1:j-1
            if x_star[i,j] < threshold
                push!(rows, i)
                push!(cols, j)
            end
        end
    end
    A_derived = sparse(rows, cols, ones(length(rows)), n, n)
    A_derived = A_derived + A_derived'
    return A_derived
end

"""
***Pulled from https://github.com/Vedangi/FastLamCC***

New type of rounding technique 

Rounds a lower bound for LambdaSTC into a feasible solution
    for a LambdaCC clustering.

This works for:
    * rounding the LambdaCC LP relaxation
    * rounding the LambdaSTC LP relaxation for lambda>=0.5 (7 - 2/lambda approx)
    * rounding the LambdaSTC LP relaxation for lambda < 0.5 ((1+lambda)/lambda - approx)

Input:
    A = adjacency matrix

    Elist = linear ordering of node-pairs

    soln = feasible solution to the CC LP relaxation
            or the LambdaSTC LP relaxation
            or the LambdaSTC ILP

Output:
    Approximation for LambdaCC.

    Step 1: Construct derived graph based on the new rounding technique
    Step 2: Apply pivot to the new graph
    Step 3: Compute LambdaCC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot 
    method on the derived graph.
"""

function CL_round_LP_new(A,Elist,soln,pivtimes,lam)
    n = size(A,1)
   
    # Elist is a linearization of all node pairs in A,
    # not just the edges in A.

    tic = time()
    tol = 1e-8
    thr = (2*lam)/((7*lam)-2)
    
    # Select the elements of soln that correspond to edges/non-edges in A according to the rounding technique
    keep = Int64[]

    if (lam < 0.5)
        thr = (lam)/(lam+1)
        for k in 1:size(Elist, 1)
            i = Elist[k,1]
            j = Elist[k,2]
            if i < j
                if ((A[i, j] != 0) || (soln[k] < (thr - tol)))  # Find the indices of the non-edge pairs that are less than thr-tol
                    push!(keep, k)
                end
            end
        end

    else
        for k in 1:size(Elist, 1)
            i = Elist[k,1]
            j = Elist[k,2]
            if i < j && soln[k] < (thr - tol)  # Find the indices of the edge pairs that are less than thr-tol
                if A[i, j] != 0
                    push!(keep, k)
                end
            end
        end
    end

   
    # These all become edges in a new derived graph Anew
    Inew = Elist[keep,1]
    Jnew = Elist[keep,2]
    Anew = sparse(Inew,Jnew,ones(length(Jnew)),n,n)
    Anew = Anew + Anew'

    Is, Js = findnz(triu(A))
    ElistA = [Is Js]
    m = sum(A)/2

    NeighbsNew = ConstructAdj(Anew,n)
    setuptime = time()-tic
    # Apply pivot on this graph multiple times,
    # returning the best output

    clus, Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
    cl_obj = check_cl_obj_faster(ElistA,clus,m,Sbin,lam)
  
    objs = cl_obj
    tic = time()
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times

        clusnew,Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
        objnew = check_cl_obj_faster(ElistA,clusnew,m,Sbin,lam)

        objs += objnew
        if objnew < cl_obj
            cl_obj = objnew
            clus = clusnew
        end
    end
    totaltime = time()-tic
    avg_time = totaltime/pivtimes + setuptime
    avg_obj = objs/pivtimes
    return round(Int64,cl_obj), clus, avg_obj, avg_time
end

"""
Fastest version of the pivot algorithm
"""
function permutation_pivot_fastest(A,Neighbs)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    clusnum = 1

    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return c, round(Int64,Sbin)
end

"""
Faster checking of LambdaCC correlation cluster objective.

Elist here is the actual edge list of a graph, don't double count.

It's not usually all pairs (unless this is a complete graph)
"""

function check_cl_obj_faster(Elist,c,m,Sbin,lam)
    n = length(c)
    pos_mis = 0
    for k = 1:size(Elist,1)
        if c[Elist[k,1]] != c[Elist[k,2]]
            pos_mis += 1
        end
    end

    neg_mis = pos_mis - m + Sbin
    return round(Int64,lam*neg_mis) + ((1-lam)*pos_mis)
end

"""
For graph with adjacency matrix A, extract triplet indices
for open wedges and triangles.

For each entry in W, the first index is the center of the wedge
each triangle is listed 3 times, once for each "center"
"""
function get_wedges_triangles(A)
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    T = Vector{Tuple{Int,Int,Int}}()
    W = Vector{Tuple{Int,Int,Int}}()
    for i = 1:n
        N = Neighbs[i]
        for jj = 1:length(N)
            j = N[jj]
            for kk = jj+1:length(N)
                k = N[kk]
                if A[j,k] == 1
                    push!(T,(i,j,k))
                else
                    push!(W,(i,j,k))
                end
            end
        end
    end
    return T, W
end