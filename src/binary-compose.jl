"""
    _compose(T₁ :: WFST{N₁, Li, L, W₁},
             T₂ :: WFST{N₂, L, Lo, W₂}) -> T :: WFST{Tuple{N₁,N₂}, Li, Lo, W}

Compose `T₁` and `T₂` by matching the output labels of `T₁` with the input labels of `T₂` lazily. 
This process is similar to matrix multiplication.

!!! note
    The details of the states of `T` are kept in `Tuple{N₁,N₂}`, this function is for demonstration purposes.

# Math
`T[x,y] = ⨁ ₖ (T₁[x,k] ⨀ T₂[k,y])` just like matrix multiplication \n
`[AB]ₘₙ = ∑ₖ Aₘₖ * Bₖₙ`
"""
function _compose(T₁::WFST{N₁, Li, L, W₁},
                  T₂::WFST{N₂, L, Lo, W₂}) where {N₁, N₂, Li, L, Lo, W₁, W₂}
    N = Tuple{N₁,N₂}
    W = promote_type(W₁, W₂)
    T = WFST{N, Li, Lo, W}()

    E₁ = T₁.states
    I₁ = T₁.starts
    F₁ = T₁.finals
    E₂ = T₂.states
    I₂ = T₂.starts
    F₂ = T₂.finals

    Q = Vector{N}()    # Vector as a Queue
    S = Set{N}()       # States Set for checking existence

    #= ONLY the initial states with composable arcs 
    are collected for the subsequent while loop =#
    for (s₁, w₁) ∈ I₁, (s₂, w₂) ∈ I₂
        for e₁ ∈ E₁[s₁], e₂ ∈ E₂[s₂]
            if ismatched(e₁, e₂)
                s₁s₂ = tuple(s₁, s₂)
                addstart(T, s₁s₂, w₁ * w₂)
                if s₁s₂ ∉ S
                    push!(Q, s₁s₂)
                    push!(S, s₁s₂)
                end
            end
        end
    end

    while !isempty(Q)
        s₁, s₂ = s₁s₂ = popfirst!(Q)
        if haskey(F₁, s₁) && haskey(F₂, s₂)
            addfinal(T, s₁s₂, F₁[s₁] * F₂[s₂])
        end

        #= states with no leaving arcs are still recorded with an empty Arc
        Vector, otherwise we would HAVE TO add the below line of code:
        `(!haskey(E₁,s₁) || !haskey(E₂,s₂)) && continue`
        here to avoid iterating into a non-exist Arc Vector in the following code =#
        for e₁ ∈ E₁[s₁], e₂ ∈ E₂[s₂]
            if olabel(e₁) == ilabel(e₂)
                n₁n₂ = tuple(next(e₁), next(e₂))
                if n₁n₂ ∉ S
                    push!(Q, n₁n₂)
                    push!(S, n₁n₂)
                end
                addarc(T, s₁s₂, n₁n₂, ilabel(e₁), olabel(e₂), weight(e₁) * weight(e₂))
            end
        end
    end # while
    return T
end



"""
    compose(T₁ :: WFST{N₁, Li, L, W₁},
            T₂ :: WFST{N₂, L, Lo, W₂}) -> T :: WFST{Int, Li, Lo, W}

Compose `T₁` and `T₂` by matching the output labels of `T₁` with the input labels of `T₂` lazily. 
This process is similar to matrix multiplication.
# Math
`T[x,y] = ⨁ ₖ (T₁[x,k] ⨀ T₂[k,y])` just like matrix multiplication \n
`[AB]ₘₙ = ∑ₖ Aₘₖ * Bₖₙ`
"""
function compose(T₁::WFST{N₁, Li, L, W₁},
                 T₂::WFST{N₂, L, Lo, W₂}) where {N₁, N₂, Li, L, Lo, W₁, W₂}
    N = Tuple{N₁,N₂}
    W = promote_type(W₁, W₂)
    T = WFST{Int, Li, Lo, W}()

    E₁ = T₁.states
    I₁ = T₁.starts
    F₁ = T₁.finals
    E₂ = T₂.states
    I₂ = T₂.starts
    F₂ = T₂.finals

    indx = -1          # state index of T
    D = Dict{N,Int}()  # Tuple{N₁,N₂} => Int
    Q = Vector{N}()    # Vector as a Queue
    S = Set{N}()       # States Set for checking existence

    #= ONLY the initial states with composable arcs 
    are collected for the subsequent while loop =#
    for (s₁, w₁) ∈ I₁, (s₂, w₂) ∈ I₂
        for e₁ ∈ E₁[s₁], e₂ ∈ E₂[s₂]
            if ismatched(e₁, e₂)
                s₁s₂ = tuple(s₁, s₂)
                if !haskey(D, s₁s₂)
                    indx = indx + 1
                    D[s₁s₂] = indx
                end
                addstart(T, D[s₁s₂], w₁ * w₂)
                if s₁s₂ ∉ S
                    push!(Q, s₁s₂)
                    push!(S, s₁s₂)
                end
            end
        end
    end

    while !isempty(Q)
        s₁, s₂ = s₁s₂ = popfirst!(Q)
        if haskey(F₁, s₁) && haskey(F₂, s₂)
            addfinal(T, D[s₁s₂], F₁[s₁] * F₂[s₂])
        end

        #= states with no leaving arcs are still recorded with an empty Arc
        Vector, otherwise we would HAVE TO add the below line of code:
        `(!haskey(E₁,s₁) || !haskey(E₂,s₂)) && continue`
        here to avoid iterating into a non-exist Arc Vector in the following code =#
        for e₁ ∈ E₁[s₁], e₂ ∈ E₂[s₂]
            if olabel(e₁) == ilabel(e₂)
                n₁n₂ = tuple(next(e₁), next(e₂))
                if !haskey(D, n₁n₂)
                    indx = indx + 1
                    D[n₁n₂] = indx
                end
                if n₁n₂ ∉ S
                    push!(Q, n₁n₂)
                    push!(S, n₁n₂)
                end
                addarc(T, D[s₁s₂], D[n₁n₂], ilabel(e₁), olabel(e₂), weight(e₁) * weight(e₂))
            end
        end
    end # while
    return T
end



"""
    ∘(T₁ :: WFST, T₂ :: WFST) -> T :: WFST
Compose two WFSTs, `T₁`'s input labels shall match `T₂`'s output labels.
"""
function Base.:∘(T₁::WFST, T₂::WFST)
    return compose(T₁, T₂)
end



