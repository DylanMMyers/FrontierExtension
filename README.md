# FrontierExtension
My implementation of the Frontier Extension algorithm, which was created by the authors of the paper "Graph Clustering in All Parameter Regimes" (https://drops.dagstuhl.de/storage/00lipics/lipics-vol170-mfcs2020/LIPIcs.MFCS.2020.39/LIPIcs.MFCS.2020.39.pdf).

Uses existing algorithms for solving LambdaCC and rounding approximated clusterings from this repo: https://github.com/nveldt/LambdaCC-Projects

# Usage
I've tried to make it as easy to run as possible, there are a few flags added to 
make it straightforward:

**--graph <_dir_>** - pass in the directory to your graph in relation to root of directory (default is the karate graph).

**--epsi <_val_>** - pass in the desired epsilon value (default is 0).

**--lfag** - determine if you want to determine the left bound of all considered lambda
values (by default and in the paper, we only compute the right bound for efficiency). Note this will increase the runtime.

**--no-prune** - Disables pruning of the set of clusterings/lambdas returned by FE.

**--wedge-constraint** - Use the open wedge constraint (default is to generate all triangle inequalities)

**--show-clusters** - Off by default as it is a lot of output, if you want to see the clusters (produce pdf visualizations of them) use this flag

# Example
**Note that the order of the commands doesn't matter**

julia main.jl --graph graphs/small/Karate.mat --epsi 1.0

julia main.jl --graph graphs/small/Karate.mat --epsi 0.5 --lflag

julia main.jl --graph graphs/medium/Karate.mat --epsi 0.1 --lflag --no-prune --wedge-constraint --show-clusters