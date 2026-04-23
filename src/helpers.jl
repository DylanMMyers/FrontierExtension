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