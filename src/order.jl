function keyorder(T::WFST{N}, bias::Int=0) where N
    d = Dict{N, Int}()
    i = 0
    for s ∈ keys(T.states)
        d[s] = i + bias
        i += 1
    end
    return d
end

"""
    reorder(T::WFST{N,I,O,W}, bias::Int=0) -> Q::WFST{Int,I,O,W}

Reorder the states of `T` by assigning new integer state indices, starting from a specified `bias`. 
The resulting WFST `Q` would have the same transitions, input/output labels, and weights as `T`. 
This is ofen used after some operations like compose, determine etc.
# Example
```julia
T = WFST((0,1),(1,2),"XXX","XXX",1)   # Create a empty WFST by args examples
addarc(T, (0,1),(10,2),"in","out", 1) # add an arc
Q = reorder(T, 10)                    # state index starting from 10
```
"""
function reorder(T::WFST{N,I,O,W}, bias::Int=0) where {N,I,O,W}
    Q = WFST{Int,I,O,W}()
    d = keyorder(T, bias)

    for (s, w) ∈ T.starts
        idx = d[s]
        addstart(Q, idx, w)
    end
    for (f, w) ∈ T.finals
        idx = d[f]
        addfinal(Q, idx, w)
    end
    for (state, arcs) ∈ T.states
        for arc ∈ arcs
            srcid = d[state]
            dstid = d[next(arc)]
            addarc(Q, srcid, dstid, ilabel(arc), olabel(arc), weight(arc))
        end
    end
    return Q
end
