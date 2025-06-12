"""
    symdists(T::WFST{N}, p::N) -> d :: Dict{N,Int}

Return `d` as the minimum consumption of symbol`s` from state `p` to all states of `T`. 
Note that the symbol-distances from `p` to unreachable states are ``∞``. See also `dists`. 
This algorithm is helpful for 
+ minimization of an DFA 
+ showing the FST in an ordered manner hierarchically. 
# Paper Ref
*A Weight Pushing Algorithm for Large Vocabulary Speech Recognition*

!!! note
    For cyclic WFST, this algorithm is also effective.

#TODO maybe pure epsilion arcs (i.e. ϵ:ϵ/w) shall be considered.
"""
function symdists(T::WFST{N}, p::N) where N
    d = Dict{N,Int}()
    r = Dict{N,Int}()
    ∞ = typemax(Int)  # very big integer pretending as infinity
    for q ∈ keys(T.states)
        d[q] = ∞      # shortest symbol distance between p and q
        r[q] = ∞      # distance added to d[q] since last time q was extracted from S
    end
    d[p] = 0          # the distance from p to p is apparently zero
    r[p] = 0          # no distance is added from p to p

    E = T.states
    S = Set{N}()
    push!(S, p)
    while !isempty(S)
        q = pop!(S) 
        R = r[q]
        r[q] = 0
        # q ──(e)──► n
        for e ∈ E[q]
            n  = next(e)
            dn = d[n]
            Rc = R + 1
            if dn ≠ min(Rc, dn)
                d[n] = min(Rc,   dn)
                r[n] = min(Rc, r[n])
                (n ∉ S) && push!(S, n)
            end
        end
    end
    return d
end


"""
    depthfwd(T::WFST{N}) -> d :: Dict{N,Int}
Return `d` as the minimum consumption of symbol`s` from start states to all states of `T`. 
"""
function depthfwd(T::WFST{N}) where N
    d = Dict{N,Int}()
    r = Dict{N,Int}()
    ∞ = typemax(Int)
    for q ∈ keys(T.states)
        d[q] = ∞  # shortest symbol distance between p and q
        r[q] = ∞  # distance added to d[q] since last time q was extracted from S
    end

    S = Set{N}()
    for start ∈ keys(T.starts)
        d[start] = 0
        r[start] = 0
        push!(S, start)
    end

    E = T.states
    while !isempty(S)
        q = pop!(S) 
        R = r[q]
        r[q] = 0
        # q ──(e)──► n
        for arc ∈ E[q]
            n  = next(arc)
            dn = d[n]
            Rc = R + 1
            if dn ≠ min(Rc, dn)
                d[n] = min(Rc,   dn)
                r[n] = min(Rc, r[n])
                (n ∉ S) && push!(S, n)
            end
        end
    end
    return d
end


"""
    depthbwd(T::WFST{N}) -> d :: Dict{N,Int}
Return `d` as the minimum consumption of symbol`s` from final states to all states of `T`. 
"""
function depthbwd(T::WFST{N}) where N
    d = Dict{N,Int}()
    r = Dict{N,Int}()
    ∞ = typemax(Int)
    for q ∈ keys(T.states)
        d[q] = ∞  # shortest symbol distance between p and q
        r[q] = ∞  # distance added to d[q] since last time q was extracted from S
    end

    S = Set{N}()
    for final ∈ keys(T.finals)
        d[final] = 0
        r[final] = 0
        push!(S, final)
    end

    #= reverse arc and just save state info
    i.e.  A ───► B     ⇒     B ────► A  =#
    E = Dict{N, Vector{N}}()
    for (src, arcs) ∈ T.states
        for e ∈ arcs
            dst = next(e)
            keyvecpush!(E, dst, src)
        end
        if !haskey(E, src)
            E[src] = Vector{N}()
        end
    end

    while !isempty(S)
        q = pop!(S) 
        R = r[q]
        r[q] = 0
        # q ──(e)──► n
        for n ∈ E[q]
            dn = d[n]
            Rc = R + 1
            if dn ≠ min(Rc, dn)
                d[n] = min(Rc,   dn)
                r[n] = min(Rc, r[n])
                (n ∉ S) && push!(S, n)
            end
        end
    end
    return d
end


"""
    depth(T::WFST{N}, direction::String) -> d :: Dict{N,Int}
Return `d` as the minimum consumption of symbol`s` from:
+ start states to all states of `T` if `direction` is "fwd"
+ final states to all states of `T` if `direction` is "bwd"
"""
function depth(T::WFST, direction::String)
    isequal(direction, "fwd") && return depthfwd(T)
    isequal(direction, "bwd") && return depthbwd(T)
    error("arg `direction` is either \"fwd\" or \"bwd\"")
end


"""
depthsets(T::WFST{N}, direction::String) -> (depth2states::Dict{Int,Vector{N}}, maxdepth::Int)
"""
function depthsets(T::WFST{N}, direction::String) where N
    MAX = 0
    keyset = Dict{Int,Set{N}}()
    for (state, dϵpth) ∈ depth(T, direction)
        keysetpush!(keyset, dϵpth, state)
        (MAX < dϵpth) && (MAX = dϵpth)
    end
    return keyset, MAX
end


"""
depthvecs(T::WFST{N}, direction::String) -> (depth2states::Dict{Int,Vector{N}}, maxdepth::Int)
"""
function depthvecs(T::WFST{N}, direction::String) where N
    MAX = 0
    keyvec = Dict{Int,Vector{N}}()
    for (state, dϵpth) ∈ depth(T, direction)
        keyvecpush!(keyvec, dϵpth, state)
        (MAX < dϵpth) && (MAX = dϵpth)
    end
    return keyvec, MAX
end

