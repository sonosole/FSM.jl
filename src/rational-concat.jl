"""
    *(T₁::WFST{N₁, I, O, W₁}, T₂::WFST{N₂, I, O, W₂}) -> T::WFST{Int, I, O, promote_type(W₁, W₂)}

Concatenate two WFSTs `T₁` and `T₂` sequentially. The resulting WFST `T` contains:
1. all arcs and states from `T₁` and `T₂`, with new index mapping to ensure no state conflicts.
2. connections from final states of `T₁` to start states of `T₂` through ϵ-labeled arcs, weighted by the product of their respective final and start state weights.
"""
function Base.:*(T₁::WFST{N₁, I, O, W₁},
                 T₂::WFST{N₂, I, O, W₂}) where {N₁, N₂, I, O, W₁, W₂}
    W = promote_type(W₁, W₂)
    T = WFST{Int, I, O, W}()

    #= states of `T₁` and `T₂` are relabeled using `keyorder` 
    to create unique state indices for each WFST. =#
    d₁ = keyorder(T₁, 0)            # idx: 0 ~ (|T₁|ₛ - 1)
    d₂ = keyorder(T₂, nstates(T₁))  # idx: |T₁|ₛ ~ (|T₁|ₛ + |T₂|ₛ - 1)

    #= make new arcs from T₁ and T₂ with new idx arcs 
    from `T₁` are added first, followed by arcs from `T₂`. =#
    for (s, arcs) ∈ T₁.states, e ∈ arcs
        addarc(T, d₁[s], d₁[next(e)], ilabel(e), olabel(e), weight(e))
    end
    for (s, arcs) ∈ T₂.states, e ∈ arcs
        addarc(T, d₂[s], d₂[next(e)], ilabel(e), olabel(e), weight(e))
    end

    # starts of T₁ → starts of T
    for (s, w) ∈ T₁.starts
        addstart(T, d₁[s], W(w))
    end
    # finals of T₂ → finals of T
    for (f, w) ∈ T₂.finals
        addfinal(T, d₂[f], W(w))
    end

    #= final states of `T₁` are connected to start states of `T₂` with 
    ϵ-transitions. so we have new arcs: {d₁[f]───ϵ:ϵ/(wf*ws)───►d₂[s]} =#
    for (f, wf) ∈ T₁.finals
        for (s, ws) ∈ T₂.starts
            addarc(T, d₁[f], d₂[s], one(I), one(O), wf * ws)
        end
    end

    return T
end


"""
    cat(T₁::WFST{N₁, I, O, W₁},
        T₂::WFST{N₂, I, O, W₂}) -> T::WFST{Int, I, O, promote_type(W₁, W₂)}

Concatenate two WFSTs `T₁` and `T₂` sequentially. The resulting WFST `T` contains:
1. all arcs and states from `T₁` and `T₂`, with new index mapping to ensure no state conflicts.
2. connections from final states of `T₁` to start states of `T₂` through ϵ-labeled arcs, weighted by the product of their respective final and start state weights.
"""
function Base.cat(T₁::WFST{N₁, I, O, W₁},
                  T₂::WFST{N₂, I, O, W₂}) where {N₁, N₂, I, O, W₁, W₂}
    return T₁ * T₂
end

