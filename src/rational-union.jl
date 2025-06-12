"""
    T::WFST{Int, I, O, W} = T₁::WFST{N₁, I, O, W₁} + T₂::WFST{N₂, I, O, W₂}

union `T₁` and `T₂` in parallel. The new WFST `T` recognizes any sequence recognized by either `T₁` or `T₂`. 
This is achieved by introducing a new start state that has ϵ-labeled arcs to the start states of both `T₁` and `T₂`. 
`W` = promote_type(`W₁`, `W₂`)
"""
function Base.:+(T₁::WFST{N₁, I, O, W₁},
                 T₂::WFST{N₂, I, O, W₂}) where {N₁, N₂, I, O, W₁, W₂}
    W = promote_type(W₁, W₂)
    T = WFST{Int, I, O, W}()

    d₁ = keyorder(T₁, 0)            # idx: 0 ~ (|T₁|ₛ - 1)
    d₂ = keyorder(T₂, nstates(T₁))  # idx: |T₁|ₛ ~ (|T₁|ₛ + |T₂|ₛ - 1)

    for (s, arcs) ∈ T₁.states
        for e ∈ arcs
            addarc(T, d₁[s], d₁[next(e)], ilabel(e), olabel(e), weight(e))
        end
    end
    for (s, arcs) ∈ T₂.states
        for e ∈ arcs
            addarc(T, d₂[s], d₂[next(e)], ilabel(e), olabel(e), weight(e))
        end
    end

    #= introduce a new start state L that has ϵ-labeled
    arcs connected to old starts of T₁ and T₂ =#
    L = nstates(T₁) + nstates(T₂); addstart(T, L)
    for (s, w) ∈ T₁.starts addarc(T, L, d₁[s], one(I), one(O), w) end
    for (s, w) ∈ T₂.starts addarc(T, L, d₂[s], one(I), one(O), w) end

    # (old finals of T₁ and T₂) become (new finals of T)
    for (f, w) ∈ T₁.finals addfinal(T, d₁[f], w) end
    for (f, w) ∈ T₂.finals addfinal(T, d₂[f], w) end
    return T
end


"""
    union(T₁::WFST{N₁, I, O, W₁},
          T₂::WFST{N₂, I, O, W₂}) -> T::WFST{Int, I, O, promote_type(W₁, W₂)}

union `T₁` and `T₂` in parallel. The new WFST `T` recognizes any sequence recognized by either `T₁` or `T₂`. 
This is achieved by introducing a new start state that has ϵ-labeled arcs to the start states of both `T₁` and `T₂`.
"""
function Base.union(T₁::WFST{N₁, I, O, W₁},
                    T₂::WFST{N₂, I, O, W₂}) where {N₁, N₂, I, O, W₁, W₂}
    return T₁ + T₂
end


"""
    T::WFST{Int, I, O, W} = T₁::WFST{N₁, I, O, W₁} ∪ T₂::WFST{N₂, I, O, W₂}

union `T₁` and `T₂` in parallel. The new WFST `T` recognizes any sequence recognized by either `T₁` or `T₂`. 
This is achieved by introducing a new start state that has ϵ-labeled arcs to the start states of both `T₁` and `T₂`. 
`W` = promote_type(`W₁`, `W₂`)
"""
function ∪(T₁ :: WFST{N₁, I, O, W₁},
            T₂ :: WFST{N₂, I, O, W₂}) where {N₁, N₂, I, O, W₁, W₂}
    return T₁ + T₂
end

