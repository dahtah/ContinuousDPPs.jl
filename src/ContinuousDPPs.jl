"""
Placeholder for a short summary about ContinuousDPPs.
"""
module ContinuousDPPs

using Determinantal,StaticArrays,LinearAlgebra,LowRankApprox,Printf

include("domains.jl")
include("quadrature.jl")
include("nystrom.jl")
include("dpp.jl")
export Domain,QuadDomain,integrate,compress,NystromApp,eigfun,ContinuousDPP,intensity
end # module
