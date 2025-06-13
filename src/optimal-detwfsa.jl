"""
    detwfsa(T::WFST{N, I, I, W}) -> D::WFST{Vector{Tuple{N,W}}, I, I, W}

Determinize `T`, producing a new deterministic **WFSA** `D`. A deterministic WFST has **ONLY** one path 
for a given input sequence. In a nondeterministic WFST, a state with multiple leaving arcs having 
the same input label introduces ambiguity. The determinization process eliminates this ambiguity 
by ensuring that each input symbol corresponds to a unique leaving arc at any state. After `detwfsa`, 
calling `reorder` on `D` is recommoded to re-map the state to integer indices. It has an alias `adet`.

!!! note
    If `T` contains ϵ:ϵ-labeled transitions, ensure they are removed before calling `detwfsa`.

# Design
    ───────────────────────────────────────────────────────────────────────
    p′ = pvs = {(p,v)}, a set of weighted source states
    q′ = qvs = {(q,v)}, a set of weighted destination states
    where    p denotes source state,
             q denotes destination state, and
             v denotes the residual weight of p or q.
    ───────────────────────────────────────────────────────────────────────
    contains two main loops:
        ● 1-st level: iteration of the weighted set of states (p,v) ∈ Q[p′]
        ● 2-nd level: e ∈ E(p)
    ───────────────────────────────────────────────────────────────────────
"""
function detwfsa(T::WFST{N, I, I, W}) where {N, I, W}
    T = onestart(T)
    Ṅ = Vector{Tuple{N,W}}
    D = WFST{Ṅ, I, I, W}()
    S = Set{Ṅ}()     # set of weighted state Set for existence checking
    Q = Vector{Ṅ}()  # queue for doing visiting, Vector struct is fine
    for (p, v) ∈ T.starts
        pv = [(p, v)]
        push!(Q, pv)
        addstart(D, pv)
    end

    E = T.states  # state => edges
    F = T.finals  # state => weight
    while !isempty(Q)
        pvs = pop!(Q)
        xw  = Dict{I, W}()
        xqw = Dict{I, Dict{N,W}}()

        for (p, v) ∈ pvs, e ∈ E[p]
            x = ilabel(e)
            w = weight(e)
            q =   next(e)
            w̌ = v * w
            # ⊕ x:x/w info
            if haskey(xw, x)
                xw[x] += w̌
            else
                xw[x] = w̌
            end
            # ⊕ x:x/w ─► q info
            if haskey(xqw, x)
                if haskey(xqw[x], q)
                    xqw[x][q] += w̌
                else
                    xqw[x][q] = w̌
                end
            else
                xqw[x] = Dict{N,W}(q => w̌)
            end
        end

        for (x, w) ∈ xw
            # qvs is the dst states' container with residual weight
            qvs = Vector{Tuple{N,W}}()
            for (q, w̌) ∈ xqw[x]
                #= residual weight w⁻¹ * w̌, where w is the `prefix` weight but w⁻¹ is not
                pre-computed, because its type may be an inexpressible type, e.g String =#
                push!(qvs, (q, w̌ / w))
            end
            addarc(D, pvs, qvs, x, x, w)

            if qvs ∉ S
                push!(S, qvs)
                push!(Q, qvs)
            end

            hasfinal = false
            finalval = zero(W)
            for (q, v) ∈ qvs
                if isfinal(T, q)
                    hasfinal = true
                    finalval += v * F[q]
                end
            end
            hasfinal && addfinal(D, qvs, finalval)
        end
    end
    return D
end


"""
    isdet(T::WFST)::Bool

Return `true` if `T` is a deterministic **WFST**. A deterministic WFST requires:
+ Unique start state.
+ Any state do not share the same input label.
"""
function isdet(T::WFST{N,I}) where {N,I}
    nstarts(T) > 1 && return false
    for arcs ∈ values(T.states)
        iset = Set{I}()
        for arc ∈ arcs
            i = ilabel(arc)
            if i ∉ iset
                push!(iset, i)
            else
                return false
            end
        end
    end
    return true
end


issequential(T::WFST) = isdet(T)
isfunctional(T::WFST) = isdet(reverse(T))

@doc """
    issequential(T::WFST)::Bool
return true if and only if `T` is deterministic with respect to the input labels of their arcs.
""" issequential


@doc """
    isfunctional(T::WFST)::Bool
return true if and only if `T` has only one output symbol sequence for any given input symbol sequence.
""" isfunctional


"""
    adet(T::WFST)
`adet` is an alias of `detwfsa`
"""
adet(T::WFST) = detwfsa(T)
