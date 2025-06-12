"""
    intersect(T₁::WFST{N₁, L, L}, T₂::WFST{N₂, L, L}; isfsa::Bool=true) -> T :: WFST{{N₁,N₂}, L, L}
Compose two WFSAs.
"""
function Base.intersect(T₁::WFST{N₁, L, L}, T₂::WFST{N₂, L, L}; isfsa::Bool=true) where {N₁, N₂, L}
    isfsa && return compose(T₁, T₂)
    return compose(ooproj(T₁), iiproj(T₂))
end


"""
    ∩(T₁::WFST{N₁, L, L}, T₂::WFST{N₂, L, L}; isfsa::Bool=true) -> T :: WFST{{N₁,N₂}, L, L}
Compose two WFSAs.
"""
function Base.:∩(T₁::WFST{N₁, L, L}, T₂::WFST{N₂, L, L}; isfsa::Bool=true) where {N₁, N₂, L}
    isfsa && return compose(T₁, T₂)
    return compose(ooproj(T₁), iiproj(T₂))
end
