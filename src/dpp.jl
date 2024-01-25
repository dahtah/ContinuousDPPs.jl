struct ContinuousDPP
    na :: NystromApp
    pr :: Vector{Float64}
end


function Base.show(io::IO,pp :: ContinuousDPP)
    str = @sprintf "Continuous DPP of rank %i. Expected size %.2f" length(pp.pr) sum(pp.pr)
    print(io,str)
end


function ContinuousDPP(na :: NystromApp;mode=:L,m=-1)
    if mode == :L #L-ensemble
        α = (m > 0 ? Determinantal.solve_sp(na.λ,m) : 1.0)
        pr = @. na.λ / (na.λ + 1/α)
    elseif mode == :K #marginal kernel
        s = sum(na.λ)
        α = m/s
        if α*na.λ[1] > 1
            α = 1/na.λ[1]
            maxm = α*s
            @info "Cannot scale to required size. Setting to max. size, $(maxm)"
        end
        @show α
        pr = α*na.λ
    else
        error("Mode should be :L or :K")
    end
    ContinuousDPP(na,pr)
end

function sample_lvg(na :: NystromApp,ind)
    wlvg=sum.(eachrow(na.F[:,ind].^2)) .* na.q.wq
    x=na.q.xq[:,Determinantal.wsample(wlvg,sum(wlvg))]
    rand_lvg_mhstep(na,ind,x)
end

function sample_lvg(na :: NystromApp,ind,n)
    wlvg=Determinantal.lvg(na.F[:,ind]) .* na.q.wq
    s = sum(wlvg)
    map(1:n) do _
        x=na.q.xq[:,Determinantal.wsample(wlvg,s)]
        rand_lvg_mhstep(na,ind,x)
    end
end

function rand_lvg_mhstep(na :: NystromApp,ind,x,τ=length(na.q.wq)^(-1/length(x)))
    l1 = sum(eigfun(na,x,ind).^2)
    xp = x+τ*randn(length(x))
    if in_domain(na.q,xp)
        l2 = sum(eigfun(na,xp,ind).^2)
        rat = l2/l1
        if rand() < rat
            return xp
        end
    end
    return x
end

function intensity(pp :: ContinuousDPP, x :: AbstractVector)
    ind = findall(pp.pr .> 0)
    u = eigfun(pp.na,x,ind)'
    (u.^2)*pp.pr[ind]
end

function Determinantal.sample(pp :: ContinuousDPP)
    n = length(pp.pr)
    indy = findall(rand(n) .< pp.pr)
    m = length(indy)
    if m > 0
        p = Int(round(5*m*log1p(m)))
        S=sample_lvg(pp.na,indy,p)
        vf = (v)-> eigfun(pp.na,v,indy)
        #        @info S[1],vf(S[1])
        reduce(hcat,sample_pdpp_cont(vf,S))
    end
end

function QuadDomain(pp :: ContinuousDPP,nrep=1)
    res=map(1:nrep) do _
        xq = sample(pp)
        ll = intensity(pp,eachcol(xq))
        (xq, 1 ./ ll)
    end
    xq = reduce(hcat,map((v)->v[1],res))
    wq=reduce(vcat,map((v)->v[2],res))
    QuadDomain(xq,wq / nrep,pp.na.q.dom,pp.na.q.vol)
end






function sample_pdpp_cont(U :: Function,x0)
    m = length(U(x0[1]))
    #Initial distribution
    Q = zeros(Float64,m,m)
    inds = Vector{Vector{eltype(x0[1])}}();
    f = zeros(m)
    max_attempts = 150*m
    i0 = 1;
    f = zeros(m)
    p0=0.0;
    for ind in 1:m
        accept = false
        itm =  x0[i0];i0+=1;
        nattempts = 1
        Qv = @view Q[:,1:(ind-1)]
        while !accept
            f = U(itm)
            p0 = sum(f.^2)
            pp = p0 - sum((Qv'*f).^2)
            if (rand() < pp/p0)
                accept = true
            else
                itm =  x0[i0];i0+=1;
                nattempts+=1
            end
            if (nattempts > max_attempts)
                error("Too many A/R attempts, $(nattempts) at index $(ind). Increase rank of Nyström approximation")
            end
        end
        push!(inds,itm)
        #Gram-Schmidt
        f -= Qv*(Qv'*f)
        Q[:,ind] = f/norm(f)
    end
    #collect(inds),Q
    collect(inds)
end
