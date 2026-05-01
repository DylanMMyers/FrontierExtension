# FrontierExtension
My implementation of the Frontier Extension algorithm, which was created by the authors of the paper "Graph Clustering in All Parameter Regimes" (https://drops.dagstuhl.de/storage/00lipics/lipics-vol170-mfcs2020/LIPIcs.MFCS.2020.39/LIPIcs.MFCS.2020.39.pdf).

Uses existing algorithms for solving LambdaCC and rounding approximated clusterings from this repo: https://github.com/nveldt/LambdaCC-Projects

# Benchmarks

I tested it on various graphs with epsi -> [0.01, 0.1] (on my local machine) for runtime/family size:

- **Karate.mat (n=34, m=78)** (triangle inequalities constraint) - 3.21s/5

- **Karate.mat (n=34, m=78)** (open wedge constraint) - 1.69s/4

- **dolphins.mat (n=62, m=159)** (triangle inequalities constraint) - 16.4s/5

- **dolphins.mat (n=62, m=159)** (open wedge constraint) - 1.84s/4

- **jazzA.mat (n=198, m=2742)** (open wedge constraint) - 1775.31s/4

- **Harvard500A.mat (n=500, m=2043)** (open wedge constraint) - 7061.28/8

# Requirements
- Julia 1.12+
- Gurobi (with valid license)
- Packages: Gurobi.jl, JuMP.jl, MAT.jl, SparseArrays, LinearAlgebra, Plots

# Usage
Note that using the open wedge constraint results in lower runtimes, but can result in gaps in the covering (this is fine in practice as we change lambda by a factor of (1 + epsi) so the approx guarantee is met). Calculating the left bound and pruning obviously will fix this but results in longer runtimes as extra computation is needed (nontrivial due to having to compute ORLP twice instead of once).

I've tried to make it as easy to run as possible, there are a few flags added to 
make it straightforward and provide some flexibility:

**--graph <_dir_>** - Pass in the directory to your graph in relation to root of directory (default is the karate graph).

**--epsi <_val_>** - Pass in the desired epsilon value (default is 1).

**--no-prune** - Disables pruning of the set of clusterings/lambdas returned by FE as well as the calculating of the left bound. For the returned intervals, they will only report an accurate right bound (they implicitly cover further to the left, however for a given lambda it will be its own leftbound). Reduces runtime substantially by reducing ORLP calls from two to one.

**--wedge-constraint** - Use the open wedge constraint (default is to generate all triangle inequalities)

**--show-clusters** - Off by default as it is a lot of output, if you want to see the clusters (produce pdf visualizations of them) use this flag

# Example
**Note that the order of the commands doesn't matter**

julia main.jl --graph graphs/small/Karate.mat --epsi 1.0

julia main.jl --graph graphs/small/dolphins.mat --epsi 0.5 --wedge-constraint

julia main.jl --graph graphs/small/jazzA.mat --epsi 0.1 --no-prune --wedge-constraint --show-clusters