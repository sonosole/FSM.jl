"""
    rdrdfsa(T::WFST, inorder::Bool=true)

Brzozowski's algorithm to minimize a WFSA. The worst-case complexity of this  
algorithm is exponential in the number of states of the input automaton. This 
holds regardless of whether the input is a NFA or a DFA.
"""
function rdrdfsa(T::WFST, inorder::Bool=true)
    R = adet(reverse(adet(reverse(T))))
    if inorder
        return reorder(R)
    else
        return R
    end
end


"""
    rdrd(T::WFST, inorder::Bool=true)

Brzozowski's algorithm to minimize a WFST. The worst-case complexity of this  
algorithm is exponential in the number of states of the input automaton. This 
holds regardless of whether the input is a NFA or a DFA. Anyway, Brzozowski's 
double reversion deterministic method can merge paths with the same prefix or 
suffix very efficiently, so it's indeed a good choice for a dictionary WFST.
"""
function rdrd(T::WFST, inorder::Bool=true)
    R = tdet(reverse(tdet(reverse(T))))
    if inorder
        return reorder(R)
    else
        return R
    end
end


@inline function keysetpush!(kset::Dict{K, Set{N}}, k::K, v::N) where {K,N}
    if !haskey(kset, k)
        push!(kset,  k=>Set{N}())
    end
    return push!(kset[k], v)
end

@inline function Base.:-(x::Set{N}, y::Set{N}) where N
    return setdiff(x, y)
end

@inline function Base.min(x::Set{N}, y::Set{N}) where N
    (length(x) < length(y)) && 
    return x;
    return y;
end

@inline function sublength(x::Set{N}, y::Set{N}) where N
    return length(x) - length(y)
end


# Error Test Dict
# ```text
# aachen aa k ax n
# aachen's aa k ax n z
# aardvark aa d v aa k
# aardvarks aa d v aa k s
# ```
"""
    min(T::WFST)
Minimization of `T` via Hopcroft's algorithm. It could be applied to any 
determinized WFST including cyclic ones. So it's usually used for WFSTs 
having complex structures (e.g. n-gram language model). If the WFST is 
acyclic, we could use a more efficient algorithm called Revuz's algorithm. 
Minimization has no **lazy** implementation.

!!! note
    Before minimization, the input WFSTs shall be weight pushed.

# Pseudocode of Main Part
```julia
P = {F, Q - F}
W = {F, Q - F}
while !isempty(W) 
│   S = pop!(W)
│   for a ∈ Σ
│   │   Xa = {x | δ(x, a) ∈ S}
│   │   for B ∈ P where B1 = (B ∩ Xa) ≠ ∅, B2 = (B - B1) ≠ ∅
│   │   │   P = (P - {B}) ∪ {B1, B2}
│   │   │   if B ∈ W
│   │   │   │   W = (W - {B}) ∪ {B1, B2}
│   │   │   else
│   │   │   │   W = W ∪ {min(B1,B2)}
│   │   │   └─
│   │   └─
│   └─
└─
```
"""
function Base.min(T::WFST{N,I,O,V}) where {N,I,O,V}
    T = onefinal(T)
    P = Vector{Set{N}}()  # Partition sets
    W = Vector{Set{N}}()  # Waiting sets
    NFS = Set{N}(s for (s, _) ∈ T.states if !isfinal(T, s))  # Non Ninal States
    FIS = Set{N}(f for (f, _) ∈ T.finals                  )  # FInal States
    push!(P, FIS, NFS)    # P = {F, Q - F}
    push!(W, FIS)         # W = {min(F, Q - F)}, |F|≤|Q-F| is met after onefinal

    R = reverse(T)
    E = R.states
    while !isempty(W)
        S = pop!(W)
        Riow = Dict{Tuple{I,O,V}, Set{N}}()
        for s ∈ S, e ∈ E[s]
            a = ilabel(e), olabel(e), weight(e)
            x = next(e)
            keysetpush!(Riow, a, x)
        end
        for Xa ∈ values(Riow), B ∈ P
            #= B being splittable by Xa means: (B1 ≠ ∅, B2 ≠ ∅),  
            which is equivalent to conditions: (B ⊈ Xa, B1 ≠ ∅) =#
            B1 = intersect(B, Xa)        # B1 = B ∩ Xa ≠ ∅
            ΔL = sublength(B, B1)        # B2 = ∅ if |B2|=0
            isempty(B1) && continue      # here if B1 = ∅, then Xa can't split B
            iszero(ΔL)  && continue      # here if B2 = ∅, then B ⊆ Xa
            # ☢ do NOT in-place split B via B1 before judging B ∈ W (be careful not B ⊆ W) ☢
            if B ∈ W                    # ✉✉✉✉✉✉✉   W = {B, ...}, P = {B, ...}, 
                B2 = setdiff!(B, B1)     # B = B1 ∪ B2 ⇒ W = {B2,...}, P = {B2,...}
                push!(P, B1)             # P = {B2,B1,...}
                push!(W, B1)             # W = {B2,B1,...}
            else                         # ✉✉✉✉✉✉✉  W = {...}, P = {B, ...} 
                B2 = setdiff!(B, B1)     # B = B1 ∪ B2 ⇒ W = {...}, P = {B2,...}
                push!(P, B1)             # P = {B2,B1,...}
                push!(W, min(B1, B2))    # W = {min(B1,B2),...}
            end
        end
    end # while

    # blockidx::Int to which the state::N belongs
    id = Dict{N, Int}()
    for (blockidx, B) ∈ enumerate(P), state ∈ B
        id[state] = blockidx
    end
    M = WFST{Int, I,O,V}()
    for (f, w) ∈ T.finals addfinal(M, id[f], w) end
    for (s, w) ∈ T.starts addstart(M, id[s], w) end
    for (state, arcs) ∈ T.states, e ∈ arcs
        addarc(M, id[state], id[next(e)], ilabel(e), olabel(e), weight(e))
    end
    return uniquearcs!(M)
end


"""
    detailedmin(T::WFST{N,I,O,V}) -> M :: WFST{Vector{N},I,O,V}
Minimization of `T`. Details of the states-partition would be kept in `Vector{N}` indexed states. 
It's intended for demonstration, not for practical useage.
"""
function detailedmin(T::WFST{N,I,O,V}) where {N,I,O,V}
    T = onefinal(T)
    P = Set{Set{N}}()  # Partition set
    W = Set{Set{N}}()  # Waiting list
    NFS = Set{N}(s for (s, _) ∈ T.states if !isfinal(T, s))  # Non Ninal States
    FIS = Set{N}(f for (f, _) ∈ T.finals                  )  # FInal States
    push!(P, FIS, NFS) # P = {F, Q - F}
    push!(W, FIS)      # W = {min(F, Q - F)}

    R = reverse(T)
    E = R.states
    while !isempty(W)
        S = pop!(W)
        Riow = Dict{Tuple{I,O,V}, Set{N}}()
        for s ∈ S, e ∈ E[s]
            a = ilabel(e), olabel(e), weight(e)
            x = next(e)
            keysetpush!(Riow, a, x)
        end
        for Xa ∈ values(Riow), B ∈ P
            B1 = B ∩ Xa; isempty(B1) && continue; # B1 = B ∩ Xa ≠ ∅
            B2 = B - B1; isempty(B2) && continue; # B2 = B - B1 ≠ ∅, not efficient but clear
            setdiff!(P, [B])           # P ← P - {B}
            push!(P, B1, B2)           # P ← P ∪ {B1, B2}
            if B ∈ W
                setdiff!(W, [B])       # W ← W - {B}
                push!(W, B1, B2)       # W ← W ∪ {B1, B2}
            else
                push!(W, min(B1, B2))  # W ← W ∪ {min(B1,B2)}
            end
        end
    end # while

    # block::Vector{N} to which the state::N belongs
    id = Dict{N, Vector{N}}()
    for B ∈ P, state ∈ B
        id[state] = collect(B)
    end
    M = WFST{Vector{N}, I,O,V}()
    for (f, w) ∈ T.finals addfinal(M, id[f], w) end
    for (s, w) ∈ T.starts addstart(M, id[s], w) end
    for (state, arcs) ∈ T.states, e ∈ arcs
        addarc(M, id[state], id[next(e)], ilabel(e), olabel(e), weight(e))
    end
    return uniquearcs!(M)
end


# maybe speed up the partitioning process
function popmin!(Q::Vector{Set{N}}) where N
    idx = 0
    len = Inf
    for i = 1:length(Q)
        setlen = length(Q[i])
        if setlen < len
            len = setlen
            idx = i
        end
    end
    return popat!(Q, idx)
end


"""
    uniquearcs!(T::WFST) -> T::WFST
Removes all duplicated redundant arcs.
"""
function uniquearcs!(T::WFST{N,I,O,W}) where {N,I,O,W}
    for arcs ∈ values(T.states)
        K = Vector{Int}()
        U = Set{Tuple{N,I,O,W}}()
        for j = 1:length(arcs)
            e = arcs[j]
            niow = next(e), e.i, e.o, e.w
            if niow ∉ U
                push!(K, j)
                push!(U, niow)
            end
        end
        keepat!(arcs, K)
    end
    return T
end


# to be corrected ...
function depthmin(T::WFST{N,I,O,V}) where {N,I,O,V}
    T = onefinal(T)
    R = reverse(T)
    E = R.states
    C = Vector{Set{N}}()  # Collecter of equivalent classes
    W =    Set{Set{N}}()  # Working sets

    depth2states, maxdepth = depthsets(R, "fwd")
    push!(C, depth2states[0])
    push!(W, depth2states[0])

    for d = 1 : maxdepth
        P = Set{Set{N}}()
        push!(P, depth2states[d])
        while !isempty(W)
            S = pop!(W)
            Riow = Dict{Tuple{I,O,V}, Set{N}}()
            for s ∈ S, e ∈ E[s]
                a = ilabel(e), olabel(e), weight(e)
                x = next(e)
                keysetpush!(Riow, a, x)
            end
            for Xa ∈ values(Riow), B ∈ P
                issubset(B, Xa) && continue; B1 = B ∩ Xa;
                isempty(B1)     && continue; B2 = B - B1;
                setdiff!(P, [B])
                push!(P, B1, B2)
                if B ∈ W
                    setdiff!(W, [B])
                    push!(P, B1, B2)
                else
                    push!(W, min(B1, B2))
                end
            end
        end
        for p ∈ P
            push!(C, p)
        end
        W = P
    end

    # blockidx::Int to which the state::N belongs
    id = Dict{N, Int}()
    for (blockidx, B) ∈ enumerate(C), state ∈ B
        id[state] = blockidx
    end
    M = WFST{Int, I,O,V}()
    for (f, w) ∈ T.finals addfinal(M, id[f], w) end
    for (s, w) ∈ T.starts addstart(M, id[s], w) end
    for (state, arcs) ∈ T.states, e ∈ arcs
        addarc(M, id[state], id[next(e)], ilabel(e), olabel(e), weight(e))
    end
    return uniquearcs!(M)
end

