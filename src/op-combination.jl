"""
    ideti(T::WFST, inorder::Bool=true) -> R::WFST
where `R == invert(detwfst(invert(T)))`, if `inorder` 
is true then states of `R` are re-mapped to integers.
This operation ensures T is functional.
"""
function ideti(T::WFST, inorder::Bool=true)
    R = invert(detwfst(invert(T)))
    if inorder
        return reorder(R)
    else
        return R
    end
end


function Base.:(==)(T₁::WFST{N,I,O,W}, T₂::WFST{K,I,O,W}) where {N,K,I,O,W}
    D₁ = min(push(tdet(rmeps(T₁))))
    D₂ = min(push(tdet(rmeps(T₂))))
    size(D₁) ≠ size(D₂) && return false

    S1 = D₁.states
    S2 = D₂.states
    E1 = Set{Tuple{I,O,W}}()
    E2 = Set{Tuple{I,O,W}}()

    for arcs ∈ values(S1), e ∈ arcs
        push!(E1, (ilabel(e),olabel(e),weight(e)))
    end
    for arcs ∈ values(S2), e ∈ arcs
        push!(E2, (ilabel(e),olabel(e),weight(e)))
    end
    if E1 ≠ E2
        return false
    else
        return true
    end
end

