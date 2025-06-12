mutable struct Arc{D,I,O,W}
    next :: D  # state where this arc goes to
    i    :: I  # input label
    o    :: O  # output label
    w    :: W  # weight
    function Arc(d::D, i::I, o::O, w::W) where {D,I,O,W}
        new{D,I,O,W}(d, i, o, w)
    end
end


function Base.show(io::IO, e::Arc)
    s,i,o,w = e.next, e.i, e.o, e.w
    IL = isone(i) ? "ϵ" : i
    OL = isone(o) ? "ϵ" : o
    return print(io, "→$s $IL:$OL/$w")
end


function Base.isequal(e₁::Arc{D,I,O,W}, e₂::Arc{D,I,O,W}) where {D,I,O,W}
    s₁,i₁,o₁,w₁ = e₁.next, e₁.i, e₁.o, e₁.w
    s₂,i₂,o₂,w₂ = e₂.next, e₂.i, e₂.o, e₂.w
    !isequal(s₁,s₂) && return false
    !isequal(i₁,i₂) && return false
    !isequal(o₁,o₂) && return false
    !isequal(w₁,w₂) && return false
    return true
end

@inline function swapio!(e::Arc{D,L,L,W}) where {D,L,W}
    e.i, e.o = e.o, e.i
    return nothing
end


"""
    iseps(arc::Arc)
Return `true` if the input and output labels of `arc` are all epsilion.
"""
@inline function iseps(arc::Arc)
    return isone(ilabel(arc)) && isone(olabel(arc))
end


"""
    isieps(arc::Arc)
Return `true` if the input label of `arc` is epsilion.
"""
@inline isieps(arc::Arc) = isone(ilabel(arc))


"""
    isoeps(arc::Arc)
Return `true` if the output label of `arc` is epsilion.
"""
@inline isoeps(arc::Arc) = isone(olabel(arc))


@inline ilabel(arc::Arc) = arc.i
@inline olabel(arc::Arc) = arc.o
@inline weight(arc::Arc) = arc.w
@inline next(arc::Arc)   = arc.next


"""
    ismatched(e₁::Arc{N₁,I,L}, e₂::Arc{N₂,L,O})
Return `true` if `olabel(e₁) == ilabel(e₂)`. This definition is intend for labels of 
common types e.g. `Integer` or `String`. But if `L` is a weird type, let's say a `Tuple` 
or a `Vector`, then import `ismatched` and override it with your own label type, then 
the composition operation may create some thing cool, see the example below.

# Example
```julia
import FSM.ismatched
@inline function ismatched(e₁::Arc{D, I, Vector{Float32}},
                           e₂::Arc{D, Vector{Float32}, O}) where {D,I,O}
    similarity = cosine_similarity_metric(olabel(e₁), ilabel(e₂))
    return similarity > 0.95 ? true : false
end
```

!!! note
    The order of `e₁` and `e₂` is important, `ismatched(e₁,e₂)` does not mean `ismatched(e₂,e₁)`
"""
@inline function ismatched(e₁::Arc{N₁,I,L},
                           e₂::Arc{N₂,L,O}) where {N₁,N₂,I,L,O}
    return isequal(e₁.o, e₂.i)
end

