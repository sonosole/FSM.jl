"""
    dists(T::WFST{N,I,O,W}, p::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) -> d::Dict{N,W}

Return `d` as the shortest distance`s` from state `p` to all states of `T`. `atol` or `rtol` is 
used to control the accuracy. Note that the distances from `p` to unreachable states are zero(`W`).

# Paper Ref
*A Weight Pushing Algorithm for Large Vocabulary Speech Recognition*

# Code Explanation
Just like the viterbi decoding process, this algo terminates when all its 
reachable states have been processed. There are two main steps.
1. Initialize: p is the starting state, so `δ₀[p] = one(W)`, and `δ₀[¬p] = zero(W)`
2. Iteration: `δₜ[n] = ⨁ₑ δₜ₋₁[q] ⨀ w[e]`, where `q ──(e)──► n`, well the lattice of `δ` is not actually built, 
    `δₜ₋₁` are recorded by `r`, and `δₜ` are recorded by `d`. 

!!! note
    For cyclic WFST, this algorithm is effective if the weight of WFST is guaranteed to be k-closed, so 
    small difference between weights shall be ignored to speed up the convergence, i.e. we assume 
        ``x = y`` if ``|x - y| < δ``  
    where ``δ`` is a small value that can be ignored.
"""
function dists(T::WFST{N,I,O,W}, p::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) where {N,I,O,W}
    @assert haskey(T.states, p) begin
        error("state $p is not in the given WFST")
    end
    l =  one(W)
    o = zero(W)
    d = Dict{N,W}()
    r = Dict{N,W}()
    for q ∈ keys(T.states)
        d[q] = o  # shortest distance between p and q
        r[q] = o  # weights added to d[q] since last time q was extracted from S
    end
    d[p] = l      # path's starting state has -
    r[p] = l      # - unit starting weight

    E = T.states
    S = Set{N}()
    push!(S, p)
    while !isempty(S)
        q = pop!(S) 
        R = r[q]
        r[q] = o
        # q ──(e)──► n
        for arc ∈ E[q]
            n = next(arc)
            w = weight(arc)
            dn = d[n]
            Rw = R * w
            if !isapprox(Rw + dn, dn; atol, rtol)
                d[n] += Rw
                r[n] += Rw
                (n ∉ S) && push!(S, n)
            end
        end
    end
    #= in the real semiring if `p` has self loop
    then d[p] may be changed with non-unit value 
    so we just have to reset d[p] to one again =#
    d[p] = l
    return d
end


"""
    dist(T::WFST{N}, p::N, q::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3)

Distance between state `p` and state `q`. `atol` or `rtol` is used to control the accuracy. We can assert that 
`dist(T, p, q)` equals `dist(reverse(T), q, p)`. 

!!! note
    For cyclic WFST, this algorithm is effective if the weight of WFST is guaranteed to be k-closed, so 
    small difference between weights shall be ignored to speed up the convergence, i.e. we assume 
        ``x = y`` if ``|x - y| < δ`` 
    where ``δ`` is a small value that can be ignored.
"""
function dist(T::WFST{N}, p::N, q::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) where N
    @assert haskey(T.states, p) && haskey(T.states, q) begin
        error("state $p or $q is not in the given WFST")
    end
    distp = dists(T, p; atol, rtol)
    return distp[q]
end


"""
    mindist(T::WFST{N,I,O,W}, p::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) -> d::Tuple{N,W}

Return `d` as the shortest distance from state `p` to all final states of `T`. `atol` or `rtol` is used to control the accuracy. 
"""
function mindist(T::WFST{N,I,O,W}, p::N; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) where {N,I,O,W}
    state2dist = dists(T, p; atol, rtol)
    return sum(values(state2dist))
end



"""
    dists(T::WFST{N,I,O,W}; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) -> d::Dict{N,W}

Return `d` as the shortest distance`s` from **start** states to all states of `T`. `atol` or `rtol` is used 
to control the accuracy. The distances from **start** states to unreachable states are zero(`W`).

!!! note
    For cyclic WFST, this algorithm is effective if the weight of WFST is guaranteed to be k-closed, so 
    small difference between weights shall be ignored to speed up the convergence, i.e. we assume 
        ``x = y`` if ``|x - y| < δ``  
    where ``δ`` is a small value that can be ignored.
"""
function dists(T::WFST{N,I,O,W}; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) where {N,I,O,W}
    l =  one(W)
    o = zero(W)
    d = Dict{N,W}()
    r = Dict{N,W}()
    for q ∈ keys(T.states)
        d[q] = o  # shortest distance between p and q
        r[q] = o  # weights added to d[q] since last time q was extracted from S
    end

    S = Set{N}()
    for start ∈ keys(T.starts)
        d[start] = l  # path's starting state has -
        r[start] = l  # - unit starting weight
        push!(S, start)
    end
    E = T.states

    while !isempty(S)
        q = pop!(S) 
        R = r[q]
        r[q] = o
        # q ──(e)──► n
        for arc ∈ E[q]
            n = next(arc)
            w = weight(arc)
            dn = d[n]
            Rw = R * w
            if !isapprox(Rw + dn, dn; atol, rtol)
                d[n] += Rw
                r[n] += Rw
                (n ∉ S) && push!(S, n)
            end
        end
    end
    
    for start ∈ keys(T.starts)
        d[start] = l
    end
    return d
end



"""
    mindist(T::WFST{N,I,O,W}; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) -> d::Tuple{N,W}

Return `d` as the shortest distance from **start** state to all final states of `T`. `atol` or `rtol` is used to control the accuracy. 
"""
function mindist(T::WFST{N,I,O,W}; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) where {N,I,O,W}
    state2dist = dists(T; atol, rtol)
    return sum(values(state2dist))
end

