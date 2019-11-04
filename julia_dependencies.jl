import Pkg

Pkg.add("ArgParse")
Pkg.add("Clustering")
Pkg.add("ClusterManagers")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Distances")
Pkg.add("Distributed")
Pkg.add("StatsBase")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("MathOptInterface")
Pkg.add("Distributions")
Pkg.add("LineSearches")
Pkg.add("NLsolve")
Pkg.add("Optim")

using Pkg: @pkg_str
pkg"precompile"