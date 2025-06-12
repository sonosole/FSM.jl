"""
    maxid(T::WFST)
Return the maximum index of states of `T`.
"""
function maxid(T::WFST{N,I,O,W}) where {N<:Integer,I,O,W}
    MAX = typemin(N)
    for k ∈ keys(T.states)
        (MAX < k) && (MAX = k)
    end
    return MAX
end


function maxid(T::WFST{Vector{Tuple{N,W}},I,O,W}) where {N<:Integer,I,O,W}
    MAX = typemin(N)
    for state ∈ keys(T.states), (k,v) ∈ state
        (MAX < k) && (MAX = k)
    end
    return MAX
end


function maxid(T::WFST{Vector{Tuple{N,O,W}},I,O,W}) where {N<:Integer,I,O,W}
    MAX = typemin(N)
    for state ∈ keys(T.states), (k,o,v) ∈ state
        (MAX < k) && (MAX = k)
    end
    return MAX
end


"""
    onestart(T::WFST) -> T

make `T` having **only** one initial state.

# Example
```julia
T = WFST()
addarc(T, "0 3 A:A/1")
addarc(T, "1 3 B:B/2")
addarc(T, "2 3 C:C/3")
addstart(T, 0, 2.0)
addstart(T, 1, 3.0)
addstart(T, 2, 4.0)
addfinal(T, 3)

Q = onestart(deepcopy(T))
```
"""
function onestart(T::WFST{Int,I,O,W}) where {I,O,W}
    if nstarts(T) == 1
        return T
    end
    ϵi = one(I)
    ϵo = one(O)
    newstart = maxid(T) + 1
    for (oldstart, v) ∈ T.starts
        addarc(T, newstart, oldstart, ϵi, ϵo, v)
        rmstart(T, oldstart)
    end
    addstart(T, newstart)
    return T
end


function onestart(T::WFST{Vector{Tuple{Int,W}},I,I,W}) where {I,W}
    if nstarts(T) == 1
        return T
    end
    ϵi = one(I)
    ϵo = one(I)
    newstart = [(maxid(T)+1, one(W))]
    for (oldstart, v) ∈ T.starts
        addarc(T, newstart, oldstart, ϵi, ϵo, v)
        rmstart(T, oldstart)
    end
    addstart(T, newstart)
    return T
end


function onestart(T::WFST{Vector{Tuple{Int,O,W}},I,O,W}) where {I,O,W}
    if nstarts(T) == 1
        return T
    end
    ϵi = one(I)
    ϵo = one(O)
    newstart = [(maxid(T)+1, ϵo, one(W))]
    for (oldstart, v) ∈ T.starts
        addarc(T, newstart, oldstart, ϵi, ϵo, v)
        rmstart(T, oldstart)
    end
    addstart(T, newstart)
    return T
end




"""
    onefinal(T::WFST) -> T

make `T` having **only** one final state.

# Example
```julia
T = WFST()
addarc(T, "0 1 A/1")
addarc(T, "0 2 B/2")
addarc(T, "0 3 C/3")
addstart(T, 0)
addfinal(T, 1, 2.0)
addfinal(T, 2, 3.0)
addfinal(T, 3, 4.0)

Q = onefinal(deepcopy(T))
```
"""
function onefinal(T::WFST{Int,I,O,W}) where {I,O,W}
    if nfinals(T) == 1
        return T
    end
    ϵi = one(I)
    ϵo = one(O)
    newfinal = maxid(T) + 1
    for (oldfinal, v) ∈ T.finals
        addarc(T, oldfinal, newfinal, ϵi, ϵo, v)
        rmfinal(T, oldfinal)
    end
    addfinal(T, newfinal)
    return T
end


function onefinal(T::WFST{Vector{Tuple{Int,W}},I,I,W}) where {I,W}
    if nfinals(T) == 1
        return T
    end
    ϵi = one(I)
    ϵo = one(I)
    newfinal = [(maxid(T)+1, one(W))]
    for (oldfinal, v) ∈ T.finals
        addarc(T, oldfinal, newfinal, ϵi, ϵo, v)
        rmfinal(T, oldfinal)
    end
    addfinal(T, newfinal)
    return T
end


function onefinal(T::WFST{Vector{Tuple{Int,O,W}},I,O,W}) where {I,O,W}
    if nfinals(T) == 1
        return T
    end
    ϵi = one(I)
    ϵo = one(O)
    newfinal = [(maxid(T)+1, ϵo, one(W))]
    for (oldfinal, v) ∈ T.finals
        addarc(T, oldfinal, newfinal, ϵi, ϵo, v)
        rmfinal(T, oldfinal)
    end
    addfinal(T, newfinal)
    return T
end


