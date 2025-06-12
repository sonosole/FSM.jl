@inline function addorset(d::Dict{K,V}, k::K, v::V) where {K,V}
    if !haskey(d, k)
        d[k] = v
    else
        d[k] += v
    end
end


@inline function mulorset(d::Dict{K,V}, k::K, v::V) where {K,V}
    if !haskey(d, k)
        d[k] = v
    else
        d[k] *= v
    end
end


"""
    epsdists(T::WFST{N,I,O,W}, src::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) -> d::Dict{N,W}

return the distance of `src` to all those states `𝑺` reachable from `src` by only `ϵ` arcs. Or to be more precise 

`𝑺 = {Q[Path(src,dst)] | Path(src,dst) = ϵ} - {src}`, where `dst` has at least one non-epsilon arc if `𝑺` is not empty.

This function has an alias of `ϵdists`. This algo is a slightly-modified version of single source shortest distance 
algorithm, see `dists` for details.

!!! note
    For cyclic WFST, this algorithm is effective if the weight of WFST is guaranteed to be k-closed, so 
    small difference between weights shall be ignored to speed up the convergence, i.e. we assume 
        ``x = y`` if ``|x - y| < δ`` 
    where ``δ`` is a small value that can be ignored. A too big ``δ`` may miss some states of the ϵ-closure.
"""
function epsdists(T::WFST{N,I,O,W}, src::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) where {N,I,O,W}
    𝟎 = zero(W)
    d = Dict{N,W}()
    r = Dict{N,W}()
    d[src] = one(W) # path's starting state has -
    r[src] = one(W) # - unit starting weight

    E = T.states
    M = Set{N}()    # marked to keep
    S = Set{N}()    # for iterating
    push!(S, src)
    while !isempty(S)
        q = pop!(S) 
        R = r[q]
        r[q] = 𝟎
        # q ──(e)──► n
        for arc ∈ E[q]
            if !iseps(arc)  # store the end state of -
                push!(M, q) # - every pure epsilon path
                continue
            end
            n = next(arc)
            w = weight(arc)
            addorset(d,n,𝟎)
            addorset(r,n,𝟎)
            dn = d[n]
            Rw = R * w
            if !isapprox(Rw + dn, dn; atol, rtol)
                d[n] += Rw
                r[n] += Rw
                (n ∉ S) && push!(S, n)
            end
            isfinal(T, n) && push!(M, n)
        end
    end
    
    delete!(d, src)
    for k ∈ keys(d)
        (k ∉ M) && delete!(d, k)
    end
    return d
end


"""
    ϵdists(T::WFST{N,I,O,W}, src::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) -> d::Dict{N,W}

This function is an alias of `epsdists`.
"""
const ϵdists = epsdists
