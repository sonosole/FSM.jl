"""
    connect(T::WFST) -> R::WFST
only keep the useful states and arcs that belongs to paths (`P`) 
from initial states (`I`) to final states (`F`), i.e.

``R(x,y) : π ∈ P(I,F), i[π] = x, o[π] = y``
"""
function connect(T::WFST{N,I,O,W}) where {N,I,O,W}
    fwdset = Set{Tuple{N,N,I,O,W}}()
    bwdset = Set{Tuple{N,N,I,O,W}}()
    fwditer = Vector{N}()
    bwditer = Vector{N}()
    fwdseen = Set{N}()  # states have been seen in forward  iterating
    bwdseen = Set{N}()  # states have been seen in backward iterating

    # T: (src) ── x:y/w ──► (dst)
    for src ∈ keys(T.starts)
        push!(fwditer, src)
        push!(fwdseen, src)
    end
    TEdge = T.states
    while !isempty(fwditer)
        src = pop!(fwditer)
        for e ∈ TEdge[src]
            dst = next(e)
            if dst ∉ fwdseen
                push!(fwdseen, dst)
                push!(fwditer, dst)
            end;push!(fwdset, (src, dst, ilabel(e), olabel(e), weight(e)))
        end
    end

    # Q: (src) ◄── x:y/w ── (dst)
    Q = reverse(T)
    for dst ∈ keys(Q.starts)
        push!(bwditer, dst)
        push!(bwdseen, dst)
    end
    QEdge = Q.states
    while !isempty(bwditer)
        dst = pop!(bwditer)
        for e ∈ QEdge[dst]
            src = next(e)
            if src ∉ bwdseen
                push!(bwdseen, src)
                push!(bwditer, src)
            end;push!(bwdset, (src, dst, ilabel(e), olabel(e), weight(e)))
        end
    end

    S = T.starts
    F = T.finals
    R = WFST{N,I,O,W}()
    for (src, dst, i, o, w) ∈ intersect!(fwdset, bwdset)
        addarc(R, src, dst, i, o, w)
        isstart(T, src) && addstart(R, src, S[src])
        isfinal(T, dst) && addfinal(R, dst, F[dst])
    end

    for start ∈ keys(R.starts)
        isfinal(T, start) && addfinal(R, start, F[start])
    end

    for final ∈ keys(R.finals)
        isstart(T, final) && addstart(R, final, S[final])
    end
    return R
end

