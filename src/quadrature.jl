import QuasiMonteCarlo,FastChebInterp
using NonNegLeastSquares,KernelFunctions
QMC= QuasiMonteCarlo
FCI=FastChebInterp
#Some quadrature rules
struct QuadDomain{D}
    xq :: Matrix{Float64}
    wq :: Vector{Float64}
    dom :: Domain{D}
    vol :: Float64
end

#a heuristic to find a length-scale parameter such that a DPP with m points
#is numerically feasible
function packingls(qd :: QuadDomain{D},m) where D
    (qd.vol/m)^(1/D)
end

#Construct a default QMC quadrature
function QuadDomain(dom :: Domain{D},n) where D
    s = QMC.sample(n, dom.lw, dom.up, QMC.HaltonSample())
    keep = findall(map(dom.in_domain,eachcol(s)))
    neff = length(keep)
    @info "$(neff) points remain for quadrature"
    vol = (neff/n)*bounding_volume(dom)
    QuadDomain(s[:,keep],fill(1/neff,neff),dom,vol)
end


Base.show(io::IO,qd :: QuadDomain) = print(io,"Quadrature rule in dim $(dim(qd)) with $(npoints(qd)) points")

dim(qd :: QuadDomain{D}) where D = D
domain(qd :: QuadDomain{D}) where D = qd.dom
in_domain(qd :: QuadDomain{D},x) where D = domain(qd).in_domain(x)
points(qd :: QuadDomain{D}) where D = qd.xq
weights(qd :: QuadDomain{D}) where D = qd.wq


function integrate(qd :: QuadDomain{D},f) where D
    dot(qd.wq,map(f,eachcol(qd.xq)))
end


function npoints(qd::QuadDomain)
    size(qd.xq,2)
end
function chebp(x,a,b,n)
    u = 2*(x - a) / (b-a)  - 1
    cos(n*acos(u))
end

function chebvdm_1d(x :: AbstractVector,a,b,deg)
    [chebp(x,a,b,n) for x in x, n in 0:deg]
end

function prodcols(V,α)
    @assert length(V) == length(α)
    u = V[1][:,α[1]+1]
    for i in 2:length(α)
        u .*= V[i][:,α[i]+1]
    end
    u
end

function tdcheb(qd :: QuadDomain{D},deg) where D
    V = [chebvdm_1d(qd.xq[i,:],qd.dom.lw[i],qd.dom.up[i],deg) for i in 1:D]
    reduce(hcat,[prodcols(V,α) for d in 0:deg for α in alpha_iterator(Val(D),d)])
end

function alpha_iterator(::Val{N}, s, t=()) where {N}
    N <= 1 && return ((s, t...),) # Iterator over a single Tuple
    Iterators.flatten(alpha_iterator(Val(N-1), s-i, (i, t...)) for i in 0:s)
end

function chebvdm(qd::QuadDomain{D},order) where D
    deg = @SVector fill(order,D)
    FCI.chebvandermonde([SVector{D,Float64}(z...) for z in eachcol(qd.xq)],SVector{D,Float64}(qd.dom.lw...),SVector{D,Float64}(qd.dom.up),Tuple(deg))
end

function compress(qd::QuadDomain{D},deg :: Integer,tol=.005) where D
    V=tdcheb(qd,deg)
    compress(qd,V,tol)
end

function compress(qd::QuadDomain{D},k :: Kernel,m :: Integer,tol=.005) where D
    V=kernelmatrix(k,ColVecs(qd.xq),ColVecs(qd.xq)[1:m])
    V=Matrix(qr(V).Q)
    compress(qd,V,tol)
end

function compress(qd::QuadDomain{D},V :: AbstractMatrix,tol=.005) where D
    μ=V'*qd.wq;
    n = npoints(qd)
    r = size(V,2)
    k = 8r;
    res = 2tol;
    theta= 1;
    
    while k <= n
        U = V[1:k,:]
        qq = qr(U)
        
        μp = vec(μ'/qq.R)
        w=nonneg_lsq(Matrix(qq.Q)',μp,alg=:nnls)
        res = norm(U'*w - μ)
        if (res < tol)
            ii=findall(w .> 0)
            return QuadDomain{D}(qd.xq[:,ii],w[ii],qd.dom,qd.vol)
        else
            k = Int(floor((1+theta)*k))
        end
    end
    error("Failed to find compression with tolerance $(tol)")
end
