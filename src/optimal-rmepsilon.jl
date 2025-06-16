"""
    until_not_eps_states(T::WFST{N,I,O,W}, p::N)

Return a set of states just surrounding the eps-closure of state `p`. e.g.

`` (p) ─ ϵ:ϵ/w ─► (x) ─ ϵ:ϵ/w ─► (y) ─ i:o/w ─► (q) ``

where ``q`` would be kept, ``x`` and ``y`` would be ignored.
"""
function until_not_eps_states(T::WFST{N,I,O,W}, src::N) where {N,I,O,W}
    E = T.states
    S = Set{N}()
    Q = Vector{N}()

    push!(Q, src)
    while !isempty(Q)
        p = pop!(Q)
        for e ∈ E[p]
            q = next(e)
            if isone(ilabel(e)) && isone(olabel(e))
                !isequal(p,q) && push!(Q, q)
            else
                q ∉ S && push!(S, q)
            end
        end
    end
    return S
end


"""
    ϵindegree(T::WFST{N}, src::N) -> qi::Dict{N,Int}
return the indegree of those states reachable from `src` via pure `ϵ` arcs.
Or to be more precise, all those states on ϵ-paths is in set \n
`{Q[Path(src,dst)] | Path(src,dst) = ϵ} - {src}`. 

# Diagram (best viewed in REPL)
      ┌────────────── ϵ ───────────────┐
      │                                ↓
    (src) ─── ϵ ───► (t1) ─── ϵ ───► (t2) ─── x ───► (q) 
                     ~~~~            ~~~~
    indegrees:        1               2
"""
function ϵindegree(T::WFST{N,I,O,W}, src::N) where {N,I,O,W}
    E  = T.states
    qi = Dict{N,Int}()  # Dict: state => indegrees
    qs = Vector{N}()    # state container for iterating

    push!(qs, src)
    while !isempty(qs)
        p = pop!(qs)    # or popfirst!
        for e ∈ E[p]
            if iseps(e)
                q = next(e)
                if !haskey(qi, q)
                    push!(qs, q)    # q was not in the queue, ensure enqueue ONLY once
                    qi[q] = 1       # maybe q has only one incoming ϵ-arc
                else
                    qi[q] += 1      # here q has multiple incoming ϵ-arcs
                end
            end
        end
    end
    return qi
end


"""
    ϵtoposort(T::WFST{N}, src::N)
return the topo-sorted of those states reachable from `src` via `ϵ` arcs.
Or to be more precise, all those sorted states on ϵ-paths is in set 
`{Q[Path(src,dst)] | Path(src,dst) = ϵ} - {src}`.
"""
function ϵtoposort(T::WFST{N,I,O,W}, src::N) where {N,I,O,W}
    E  = T.states
    qi = ϵindegree(T, src)  # Dict: state => indegrees
    qs = Vector{N}()        # state container for iterating
    rank = Vector{N}()      # topological ordered states

    push!(qs, src)
    while !isempty(qs)
        p = pop!(qs)
        for e ∈ E[p]
            if iseps(e)
                q = next(e)
                qi[q] -= 1
                if iszero(qi[q])
                    push!(qs, q)
                    push!(rank, q) # pushed here because src is not included in `rank`
                end
            end
        end
    end
    return rank
end


function epsqv(T::WFST{N,I,O,W}, src::N) where {N,I,O,W}
    E  = T.states
    qw̌ = Dict{N,W}()   # for returning
    qv = Dict{N,W}()   # for accumulating
    qs = Vector{N}()   # state container for iterating
    qi = ϵindegree(T, src)
    push!(qv, src=>one(W))
    push!(qs, src)
    while !isempty(qs)
        p = pop!(qs)
        v = qv[p]
        for e ∈ E[p]
            if iseps(e)
                q = next(e)
                w̌ = v * weight(e)
                # accumulating weights
                if haskey(qv, q)
                    qv[q] += w̌
                else
                    qv[q] = w̌
                end
                # and iterating in topological order
                qi[q] -= 1
                iszero(qi[q]) && push!(qs, q)
                isfinal(T, q) && push!(qw̌, q=>qv[q])
            else
                if !isequal(p, src)
                    #= arcs starting from src must be filted out! i.e. (src)── x:y/w ──►(dst), where x or y is not ϵ. so when p==src, do 
                    NOT accumulate `src=>v` into qw̌ , otherwise (src)── x:y/𝟙 ──►(src) arc would appear after ϵ-removal operation =#
                    qw̌[p] = v
                end
            end
        end
    end
    return qw̌
end


# ONLY deal with fst having no cyclic
function rm_acyclic_eps(T::WFST{N,I,O,W}) where {N,I,O,W}
    Q = WFST{N,I,O,W}()
    for (start, w) ∈ T.starts
        addstart(Q, start, w)
    end
    for (final, w) ∈ T.finals
        addfinal(Q, final, w)
    end
    𝟘  = zero(W)
    Qᶠ = Q.finals
    F  = T.finals
    E  = T.states
    for (p, arcs) ∈ T.states
        # add non ϵ:ϵ arcs
        for arc ∈ arcs
            !iseps(arc) && addarc(Q, p, arc)
        end
        for (q, w̌) ∈ epsqv(T, p)
            # ↓ discard useless path
            iszero(w̌) && continue
            # ↓ add (p)→(r) arcs
            for e ∈ E[q]
                r =   next(e)
                i = ilabel(e)
                o = olabel(e)
                w = weight(e)
                if !(isone(i) && isone(o))
                    # ϵ:ϵ arcs must be excluded
                    addarc(Q, p, r, i, o, w̌ * w)
                end
            end

            if isfinal(T, q)
                if !isfinal(Q, p)    # note: not !isfinal(T, p)
                    #= if q is a final state of T, then whether 
                    p is a final state of Q or not, p would be set as a final
                    state with zero weight for accumulating all ϵ:ϵ/w̌ paths =#
                    addfinal(Q, p, 𝟘)
                end
                Qᶠ[p] += w̌ * F[q]
            end
        end
    end
    return Q
end


# https://cs.nyu.edu/~mohri/pub/ Weighted Automata Algorithms
"""
    rmeps(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2)
Apply epsilon removal to `T`. Keyword argument `atol` or `rtol` is used to control the accuracy, 
this is especially useful when dealing with cyclic WFST. `rmϵ` is an alias of `rmeps`.

!!! note
    For cyclic WFST, this algorithm is effective if the weight of WFST is guaranteed to be k-closed, so 
    small difference between weights shall be ignored to speed up the convergence, i.e. we assume 
        ``x = y`` if ``|x - y| < δ`` 
    where ``δ`` is a small value that can be ignored. A too big ``δ`` may miss some states of the ϵ-closure.
"""
function rmeps(T::WFST{N,I,O,W}; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) where {N,I,O,W}
    Q = WFST{N,I,O,W}()
    for (start, w) ∈ T.starts
        addstart(Q, start, w)
    end
    for (final, w) ∈ T.finals
        addfinal(Q, final, w)
    end
    𝟘  = zero(W)
    Qᶠ = Q.finals
    F  = T.finals
    E  = T.states
    for (p, arcs) ∈ T.states
        # add non ϵ:ϵ arcs
        for arc ∈ arcs
            !iseps(arc) && addarc(Q, p, arc)
        end
        for (q, w̌) ∈ ϵdists(T, p; atol, rtol)
            # ↓ discard useless path
            iszero(w̌) && continue
            # ↓ repalce (p)── ⋯ ──►(q,w̌)──i:o/w──►(r) with (p)──i:o/w̌*w──►(r)
            for e ∈ E[q]
                if !iseps(e) # ϵ:ϵ arcs must be excluded
                    addarc(Q, p, next(e), ilabel(e), olabel(e), w̌ * weight(e))
                end
            end
            if isfinal(T, q)
                if !isfinal(Q, p)    # note: not !isfinal(T, p)
                    #= if q is a final state of T, then whether 
                    p is a final state of Q or not, p would be set as a final
                    state with zero weight for accumulating all ϵ:ϵ/w̌ paths =#
                    addfinal(Q, p, 𝟘)
                end
                Qᶠ[p] += w̌ * F[q]
            end
        end
    end
    return Q
end


"""
    rmϵ(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2)
alias of `rmeps`.
"""
const rmϵ = rmeps



# TODO
"""
    iepsnorm(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2)
Returns an equivalent WFST that is input epsilon-normalized, so as to ensure 
on each path all non-epsilon input labels is ahead of any epsilon input label. 
Output epsilon-normalized WFST is defined similarly. The input WFST needs to be functional.
Keyword argument `atol` or `rtol` is used to control the accuracy, this is 
especially useful when dealing with cyclic WFST. 

!!! note
    For cyclic WFST, this algorithm is effective if the weight of WFST is guaranteed to be k-closed, so 
    small difference between weights shall be ignored to speed up the convergence, i.e. we assume 
        ``x = y`` if ``|x - y| < δ`` 
    where ``δ`` is a small value that can be ignored. A too big ``δ`` may miss some states of the ϵ-closure.
    
"""
function iepsnorm(T::WFST{N,I,O,W}; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) where {N,I,O,W}
end
