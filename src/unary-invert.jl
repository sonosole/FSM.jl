"""
    swapio(T::WFST) -> Q::WFST

Swap the input and output labels of the given WFST `T`. Creates a new WFST `Q` where:
1. The input labels of each arc in `T` become the output labels in `Q`.
2. The output labels of each arc in `T` become the input labels in `Q`.
3. The start and final states of `T` remain unchanged in `Q`.
It has an alias `invert`
# Example
```julia
T = WFST()                     # create a WFST
addarc(T, 0,1,"hello","kitty") # add an arc (0) ─── "hello":"kitty"/1 ───► (1)
Q = swapio(T)                  # invert the WFST
```
"""
function swapio(T::WFST{N, I, O, W}) where {N, I, O, W}
    Q = WFST{N, I, O, W}()

    for (s, w) ∈ T.starts addstart(Q, s, w) end
    for (f, w) ∈ T.finals addfinal(Q, f, w) end

    for (src, arcs) ∈ T.states
        for arc ∈ arcs
            addarc(Q, src, next(arc), olabel(arc), ilabel(arc), weight(arc))
        end
    end
    return Q
end


"""
    swapio!(T::WFST) -> T
Inplace version of `swapio`. It has an alias `invert!`
"""
function swapio!(T::WFST{N, L, L, W}) where {N, L, W}
    for arcs ∈ values(T.states)
        for arc ∈ arcs
            swapio!(arc)
        end
    end
    return T
end


# alias
"""
    invert(T::WFST) -> Q::WFST

Swap the input and output labels of the given WFST `T`. Creates a new WFST `Q` where:
1. The input labels of each arc in `T` become the output labels in `Q`.
2. The output labels of each arc in `T` become the input labels in `Q`.
3. The start and final states of `T` remain unchanged in `Q`.
It has an alias `swapio`.

# Example
```julia
T = WFST()                     # create a WFST
addarc(T, 0,1,"hello","kitty") # add an arc (0) ─── "hello":"kitty"/1 ───► (1)
Q = invert(T)                  # invert the WFST
```
"""
invert(T::WFST{N, I, O, W}) where {N, I, O, W} = swapio(T)


"""
    invert!(T::WFST) -> T
Inplace version of `invert`. It has an alias `swapio!`
"""
invert!(T::WFST{N, L, L, W}) where {N, L, W} = swapio!(T)
