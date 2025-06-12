"""
    iiproj(T::WFST{N, I, O}) -> Q::WFST{N, I, I}
converting a transducer into an acceptor by omitting output labels
"""
function iiproj(T::WFST{N,I,O,W}) where {N,I,O,W}
    Q = WFST{N, I, I, W}()

    for (s, w) ∈ T.starts addstart(Q, s, w) end
    for (f, w) ∈ T.finals addfinal(Q, f, w) end

    for (src, arcs) ∈ T.states, arc ∈ arcs
        addarc(Q, src, next(arc), ilabel(arc), ilabel(arc), weight(arc))
    end
    return Q
end


"""
    ooproj(T::WFST{N,I,O}) -> Q::WFST{N,O,O}

converting a transducer into an acceptor by omitting input labels
"""
function ooproj(T::WFST{N,I,O,W}) where {N,I,O,W}
    Q = WFST{N, O, O, W}()

    for (s, w) ∈ T.starts addstart(Q, s, w) end
    for (f, w) ∈ T.finals addfinal(Q, f, w) end

    for (src, arcs) ∈ T.states, arc ∈ arcs
        addarc(Q, src, next(arc), olabel(arc), olabel(arc), weight(arc))
    end
    return Q
end


"""
    iiproj!(T :: WFST{N, I, O}) :: WFST{N, I, I}
inplace converting a transducer into an acceptor by omitting output labels
"""
function iiproj!(T::WFST{N, L, L, W}) where {N, L, W}
    for arcs ∈ values(T.states)
        for arc ∈ arcs
            arc.o = arc.i
        end
    end
    return T
end


"""
    ooproj!(T :: WFST{N, I, O}) :: WFST{N, O, O}

inplace converting a transducer into an acceptor by omitting input labels
"""
function ooproj!(T::WFST{N, L, L, W}) where {N, L, W}
    for arcs ∈ values(T.states)
        for arc ∈ arcs
            arc.i = arc.o
        end
    end
    return T
end
