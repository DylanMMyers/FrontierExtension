# FrontierExtension
My implementation of the Frontier Extension algorithm, which was created by the authors of the paper "Graph Clustering in All Parameter Regimes" (https://drops.dagstuhl.de/storage/00lipics/lipics-vol170-mfcs2020/LIPIcs.MFCS.2020.39/LIPIcs.MFCS.2020.39.pdf).

Uses existing algorithms for solving LambdaCC and rounding approximated clusterings from this repo: https://github.com/nveldt/LambdaCC-Projects

# Benchmarks

I tested it on various graphs with epsi=0.1 (on my local machine) for time:

- **Karate.mat (n=34, m=78)** (triangle inequalities constraint) - 3.21s/

- **Karate.mat (n=34, m=78)** (open wedge constraint) - 

- **dolphins.mat (n=62, m=159)** (triangle inequalities constraint) - 16.4s

- **dolphins.mat (n=62, m=159)** (open wedge constraint) - 1.84s

- **jazzA.mat (n=198, m=2742)** (triangle inequalities constraint) - TOOK TOO LONG (could be worth testing on grace etc.)

- **jazzA.mat (n=198, m=2742)** (open wedge constraint) - 1775.31s

# Requirements
- Julia 1.12+
- Gurobi (with valid license)
- Packages: Gurobi.jl, JuMP.jl, MAT.jl, SparseArrays, LinearAlgebra, Plots

# Usage
Note that lower epsilon values and using the wedge constraint will result in larger raw covers (important to keep in mind if you plan on not using lflag/pruning). Obviosuly, pruning will fix this but results in longer runtimes as extra computation is needed (nontrivial due to having to compute ORLP twice instead of once).

I've tried to make it as easy to run as possible, there are a few flags added to 
make it straightforward and provide some flexibility:

**--graph <_dir_>** - Pass in the directory to your graph in relation to root of directory (default is the karate graph).

**--epsi <_val_>** - Pass in the desired epsilon value (default is 1).

**--lflag** - Determine if you want to compute the left bound of all considered lambda
values. By default and in the paper we only compute the right bound for efficiency. This might result in an interval like [0.3, 0.6] while it really covers [0.1, 0.6] (obviously causing the reported bounds to be inaccurate). Note this will increase the runtime. I'm writing this like it's not Dr. Veldt reading it but in case you forgot how your own algorithm operates.

**--no-prune** - Disables pruning of the set of clusterings/lambdas returned by FE.

**--wedge-constraint** - Use the open wedge constraint (default is to generate all triangle inequalities)

**--show-clusters** - Off by default as it is a lot of output, if you want to see the clusters (produce pdf visualizations of them) use this flag

# Example
**Note that the order of the commands doesn't matter**

julia main.jl --graph graphs/small/Karate.mat --epsi 1.0

julia main.jl --graph graphs/small/Karate.mat --epsi 0.5 --lflag

julia main.jl --graph graphs/medium/Karate.mat --epsi 0.1 --lflag --no-prune --wedge-constraint --show-clusters