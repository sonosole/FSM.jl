"""
    star(T::WFST{Int}) -> Q::WFST{Int}
star or called closure, is defined as 
    `` star(T) = T⁰ + T¹ + T² + ... + Tᵏ 
               = ϵ + T¹ + T² + ... + Tᵏ ``, where ``k = ∞``
"""
function star(T::WFST{Int, I, O, W}) where {I, O, W}
    newstart = maxid(T) + 1
    Q = WFST{Int, I, O, W}()
    # 1. copy all arcs first
    for (src, arcs) ∈ T.states, e ∈ arcs
        addarc(Q, src, next(e), ilabel(e), olabel(e), weight(e))
    end
    # 2. then connect all finals to all starts.
    for start ∈ keys(T.starts), (final, wf) ∈ T.finals
        addarc(Q, final, start, one(I), one(O), wf)
    end
    # 3. finals remain the same
    for (final, w) ∈ T.finals
        addfinal(Q, final, w)
    end
    # 4. connect new start to all old starts
    for (oldstart, w) ∈ T.starts
        addarc(Q, newstart, oldstart, one(I), one(O), w)
    end
    # 5. finally the new start is also set as final because we need T⁰
    addstart(Q, newstart)
    addfinal(Q, newstart)
    return Q
end


"""
    star!(T::WFST{Int}, new_init_idx::Int=-1)
in place version of `star`
"""
function star!(T::WFST{Int, I, O}, newstart::Int=-1) where {I, O}
    #= 1. connect all finals to all starts
    more finals in practice, so stay inner =#
    for (start, ws) ∈ T.starts
        for (final, wf) ∈ T.finals
            addarc(T, final, start, one(I), one(O), wf)
        end
    end
    # 2. connect new start to all old starts
    for (oldstart, w) ∈ T.starts
        addarc(T, newstart, oldstart, one(I), one(O), w)
    end
    # 3. delete all old starts
    for oldstart ∈ keys(T.starts)
        delete!(T.starts, oldstart)
    end
    # 4. finally the new start is also set as final because we need T⁰
    addstart(T, newstart)
    addfinal(T, newstart)
    return T
end


"""
    plus(T::WFST)
plus is defined as
    `` plus(T) = T¹ + T² + T³ + ... + Tᵏ ``, where  ``k = ∞``
"""
function plus(T::WFST{N, I, O, W}) where {N, I, O, W}
    Q = WFST{N, I, O, W}()
    # 1. copy all arcs infos
    for (src, arcs) ∈ T.states
        for e ∈ arcs
            addarc(Q, src, next(e), ilabel(e), olabel(e), weight(e))
        end
    end
    # 2. then connect all finals to all starts
    for (start, ws) ∈ T.starts
        for (final, wf) ∈ T.finals
            addarc(Q, final, start, one(I), one(O), wf)
        end
        addstart(Q, start, ws)
    end
    # 2. copy finals infos
    for (final, w) ∈ T.finals
        addfinal(Q, final, w)
    end
    return Q
end


"""
    plus!(T::WFST{Int})
in place version of `plus`
"""
function plus!(T::WFST{N, I, O}) where {N, I, O}
    # just need to connect all finals to all starts
    for (final, wf) ∈ T.finals
        for start ∈ keys(T.starts)
            addarc(T, final, start, one(I), one(O), wf)
        end
    end
    return T
end
