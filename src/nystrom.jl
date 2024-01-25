using KernelFunctions

struct NystromApp{K,D}
    q :: QuadDomain{D}
    kfun :: K
    H :: Matrix{Float64}
    F :: Matrix{Float64}
    λ :: Vector{Float64}
end


function Base.show(io::IO,na :: NystromApp)
    str = @sprintf "Nyström approximation of rank %i. Eigenvalue range %e,%e. Trace %2.f. " length(na.λ) minimum(na.λ) maximum(na.λ) sum(na.λ)
    print(io,str)
end


function NystromApp(q :: QuadDomain,kfun :: Kernel,maxrank=Inf)
    C=kernelmatrix(kfun,q.xq,q.xq)
    W = Diagonal(sqrt.(q.wq))
    H = W*C*W'
    if maxrank == Inf
        eg = eigen(Symmetric(H),sortby=-)
        F =W\eg.vectors
        λ = max.(eg.values,0.0)
    else #Use LowRankApprox
        ph=pheigfact(H,rank=maxrank)
        ord = sortperm(ph.values,rev=true)
        λ = max.(ph.values[ord],0.0)
        F=W\ph.vectors[:,ord]
    end
    NystromApp(q,kfun,H,F,λ)
end

function eigfun(na::NystromApp,x :: Vector{Float64},ind)
    c = [na.kfun(x,na.q.xq[:,j])*na.q.wq[j] for j in eachindex(na.q.wq)]
    #c = kernelmatrix(na.kfun,ColVecs(na.q.xq),[x])
    #vec((na.F[:,ind]'*c) ./ na.λ[ind])
    (na.F[:,ind]'*c) ./ na.λ[ind]
end

function eigfun(na::NystromApp,x :: AbstractVector,ind)
    #c = [na.kfun(x,na.q.xq[:,j])*na.q.wq[j] for j in eachindex(na.q.wq)]
    C = Diagonal(na.q.wq)*kernelmatrix(na.kfun,ColVecs(na.q.xq),x)#*
    (na.F[:,ind]'*C) ./ na.λ[ind]
end


function check_orth(na::NystromApp{K,D},q2 :: QuadDomain,ind) where K where D
    Z=eigfun(na,ColVecs(q2.xq),ind)
    Z*Diagonal(q2.wq)*Z'
end

function lvg(na::NystromApp,x :: AbstractVector,ind)
    Determinantal.lvg(eigfun(na,x,ind)')
end

function lvg(na::NystromApp,x :: Vector{Float64},ind)
    norm(eigfun(na,x,ind))^2
end

function kdiagapp(na :: NystromApp,x :: AbstractVector)
    ind = findall(na.λ .!= 0)
    u = eigfun(na,x,ind)'
    (u.^2)*na.λ[ind]
end

function kapp(na :: NystromApp,x :: AbstractVector,y :: AbstractVector)
    u = eigfun(na,x,1:length(na.λ))
    v = eigfun(na,y,1:length(na.λ))
    u'*Diagonal(na.λ)*v
end
