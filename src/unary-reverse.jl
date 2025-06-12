"""
    reverse(T::WFST) -> Q::WFST

Creates a new WFST `Q` where:
1. All arcs that were in `T` are reversed in `Q`, meaning that each arc's direction is flipped.
2. The start states of `T` become the final states of `Q`.
3. The final states of `T` become the start states of `Q`.

# Example
```julia
T = WFST()                     # create a WFST
addarc(T, 0,1,"hello","kitty") # add an arc
Q = reverse(T)                 # reverse the WFST
```
"""
function Base.reverse(T::WFST{N, I, O, W}) where {N, I, O, W}
    Q = WFST{N, I, O, W}()
    # flip arcs' direction
    for (src, arcs) ∈ T.states, e ∈ arcs
        addarc(Q, next(e), src, ilabel(e), olabel(e), weight(e))
    end
    # starts of T become finals of Q
    for (start, w) ∈ T.starts
        addfinal(Q, start, w)
    end
    # finals of T become starts of Q
    for (final, w) ∈ T.finals
        addstart(Q, final, w)
    end
    return Q
end
