
"""
    arcsort!(T::WFST, by::Function=weight, rev::Bool=true)

Sort the arcs in an `T` per state. Arcs are first transformed with the function 
`by` which has 4 options: `ilabel`, `olabel`, `weight` or `next`. The sort order 
is determined by the `rev` kwarg. `arcsort!` is a **lazy** implementation.

# Complexity
Time : ``Q * E*log(E)`` \n
Space: ``O(E)``         \n
where ``Q`` is the number of states visited and 
      ``E`` is the number of arcs visited
"""
function arcsort!(T::WFST; by::Function=weight, rev::Bool=true)
    for arcs ∈ values(T.states)
        sort!(arcs; by, rev)
    end
end


