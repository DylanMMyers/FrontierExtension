using LinearAlgebra
using SparseArrays
using KrylovKit

"""
Given an adjacency matrix and a 2-D layout, plot the graph
"""
function display_graph(A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64},grayscale = 0.0,ms = 6,lw = 1)
    f = plot(legend=false, axis = false,grid = false,xticks = false,yticks = false)
    ei,ej,w = findnz(triu(A))
    scatter!(f,xy[:,1],xy[:,2],color = RGB(grayscale,grayscale,grayscale),markersize = ms, markerstrokecolor =  RGB(grayscale,grayscale,grayscale))
    lx = [xy[ei,1]';xy[ej,1]']
    ly = [xy[ei,2]';xy[ej,2]']
    for i = 1:length(w)
        # draws line from the first point, to the second point
        plot!(f,lx[:,i],ly[:,i],color = RGB(grayscale,grayscale,grayscale), linewidth = lw)
    end
    return f
end

"""
Given an adjacency matrix and a 2-D layout, plot the graph.
Display the edges based on Xij values between 0 and 1
"""
function display_graph_X(A::SparseMatrixCSC{Float64,Int64},X,xy::Matrix{Float64},grayscale = 0.0,ms = 6,lw = 1)
    f = plot(legend=false, axis = false,grid = false,xticks = false,yticks = false)
    f = display_graph(A,xy,grayscale,ms,lw)
    ei,ej,w = findnz(triu(X))
    scatter!(f,xy[:,1],xy[:,2],color = RGB(grayscale,grayscale,grayscale),markersize = ms, markerstrokecolor =  RGB(grayscale,grayscale,grayscale))
    lx = [xy[ei,1]';xy[ej,1]']
    ly = [xy[ei,2]';xy[ej,2]']
    for i = 1:length(w)
        # draws line from the first point, to the second point
        ii = ei[i]
        jj = ej[i]
        plot!(f,lx[:,i],ly[:,i],color = :blue, linealpha = X[ii,jj], linewidth = lw*2)
    end
    return f
end

"""
Given an adjacency matrix and a 2-D layout, plot the graph.
Display the edges based on Xij values between 0 and 1
"""
function display_graph_clustering(A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64},c,grayscale = 0.0,ms = 6,lw = 1)
    f = plot(legend=false, axis = false,grid = false,xticks = false,yticks = false)
    ei,ej,w = findnz(triu(A))
    scatter!(f,xy[:,1],xy[:,2],color = RGB(grayscale,grayscale,grayscale),markersize = ms, markerstrokecolor =  RGB(grayscale,grayscale,grayscale))
    lx = [xy[ei,1]';xy[ej,1]']
    ly = [xy[ei,2]';xy[ej,2]']
    for i = 1:length(w)
        # draws line from the first point, to the second point
        ii = ei[i]
        jj = ej[i]
        plot!(f,lx[:,i],ly[:,i],color = RGB(grayscale,grayscale,grayscale), linewidth = lw)
    end
    for i = 1:maximum(c)
        C = findall(x->x==i,c)
        scatter!(f,xy[C,1],xy[C,2],markersize = 15)
    end
    return f
end

"""
Return the Laplacian matrix for a graph with adjacency matrix A
"""
function Laplacian(A::SparseMatrixCSC)
    n = size(A,1)
    d = vec(sum(A,dims = 1))    # degree vector, make it a type "Vector"!
    D = Diagonal(d)             # turn it into degree matrix
    L = D-A                     # Laplacian matrix
    return L
end

"""
Return the normalized Laplacian matrix for a graph with adjacency matrix A
"""
function nLaplacian(A::SparseMatrixCSC,dhalf = false)
    d = vec(sum(A,dims = 2))
    Dhalf = Diagonal(d.^(-1/2))
    L = I - Dhalf*A*Dhalf
    if dhalf
        return L, Dhalf
    else
        return L
    end
end

"""
Gets xy coordinates from a graph by using a simple spectral visualization

A is an n x n adjacency matrix for a graph in sparse format

    If eigenflag = true, then this uses the function eigen, for dense graphs
"""
function spectral_xy(A, eigenflag = false)
    n = size(A,1)
    L = nLaplacian(A)
    if n < 100 || eigenflag
        Lam, V = eigen(Matrix(L))
        xy = V[:,2:3]
    else
        Vl,Vc,convinfo = eigsolve(L + n*LinearAlgebra.I, 3, :SR; tol = 1e-8, maxiter = 1000, verbosity = 0)
        xy = [Array(Vc[2]) Array(Vc[3])]
    end
    return xy
end