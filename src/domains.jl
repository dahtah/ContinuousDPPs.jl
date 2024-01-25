struct Domain{D} 
    in_domain :: Function
    lw :: Vector{Float64} #lower limit of bounding region
    up :: Vector{Float64}
end

function Domain{D}(in_domain::Function) where D
    Domain{D}(in_domain,zeros(D),ones(D))
end

function dim(dm :: Domain{D}) where D
    return D
end

Base.show(io::IO,dm :: Domain) = print(io,"Domain in dim $(dim(dm))")



function bounding_volume(dm :: Domain{D}) where D
    return prod(dm.up-dm.lw)
end

#simple A/R sampler
function Base.rand(dm :: Domain{D}; maxit=10000) where D
    accept=false
    nit=0
    dl = dm.up - dm.lw
    x = dm.lw + dl.*rand(D)
    while !accept && nit < maxit
        
        if (dm.in_domain(x))
            accept=true
        else
            x = dm.lw + dl.*rand(D)
            nit+=1
        end
    end
    if nit == maxit
        error("Too many rejections: make sure the domain is not too small relative to the bounding box")
    else
        return x
    end
end
