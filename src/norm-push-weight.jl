"""
    potential(FST::WFST{N,I,O,W}; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) where {N,I,O,W}

`W` shall be a k-closed semiring, since it is impossible to deal with infinite paths. 
The tropical semiring is 0-closed. The log semiring or prob semiring is k-closed if we 
disregard small differences between weights, in such cases `atol` or `rtol` is used to 
control the accuracy. In addition, the algorithm admits any queue discipline, a priority 
queue in which the state with the largest potential is popped first is often chosen, 
because potentials become larger and fixed earlier. 

# Math
``V[s] = ⨁ₚ w[p]``, where ``p ∈ Path(F → s)``
"""
function potential(FST::WFST{N,I,O,W}; atol::Real=0,
                                       rtol::Real= atol>0 ? 0 : 1e-2) where {N,I,O,W}
    o = zero(W)
    V = Dict{N,W}()
    r = Dict{N,W}()
    for q ∈ keys(FST.states)
        V[q] = o  # shortest distance between p and q
        r[q] = o  # weights added to d[q] since last time q was extracted from S
    end

    S = Set{N}()
    for (q, w) ∈ FST.finals
        V[q] = w
        r[q] = w
        push!(S, q)
    end

    #= reverse arc and just save state and weight info
    i.e.  A ───(i:o/w)───► B   ⇒   B ───(w)───► A  =#
    E = Dict{N, Vector{Tuple{W,N}}}()
    for (src, arcs) ∈ FST.states
        for e ∈ arcs
            dst = next(e)
            keyvecpush!(E, dst, (weight(e), src))
        end
        if !haskey(E, src)
            E[src] = Vector{Tuple{W,N}}()
        end
    end

    while !isempty(S)
        q = pop!(S) 
        R = r[q]
        r[q] = o
        for (w, n) ∈ E[q]
            Vn = V[n]
            Rw = R * w
            if !isapprox(Rw + Vn, Vn; atol, rtol)
                V[n] += Rw
                r[n] += Rw
                n ∉ S && push!(S, n)
            end
        end
    end
    return V
end



"""
    pushweight(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) -> T

Push weights towards initial states. This operation reduces searching time, 
since by looking ahead the weight, unpromissing paths can be ignored in the 
early stage of the search. `atol` or `rtol` is used to control the accuracy 
for calculating the potential.

In addition, the algorithm assumes that WFST `T` is trimmed and the weight 
semiring is weakly left-divisible and zero-sum free such as tropical and 
log semirings.
"""
function pushweight(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2)
    V = potential(T; atol, rtol)
    E = T.states
    I = T.starts
    F = T.finals

    for (q, λ) ∈ I
        I[q] = λ * V[q]
    end

    for (q, w) ∈ F
        F[q] = w / V[q]
    end

    for q ∈ keys(T.states), e ∈ E[q]
        e.w = weight(e) * V[next(e)] / V[q]
    end
    return T
end


"""
    reweight(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2) -> T

Push weights towards initial states. This operation reduces searching time, 
since by looking ahead the weight, unpromissing paths can be ignored in the 
early stage of the search. `atol` or `rtol` is used to control the accuracy 
for calculating the potential.

In addition, the algorithm assumes that WFST `T` is trimmed and the weight 
semiring is weakly left-divisible and zero-sum free such as tropical and 
log semirings.
"""
function reweight(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-2)
    return pushweight(T; atol, rtol)
end



"""
    push(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3) -> T

Push weights towards initial states. This operation reduces searching time, 
since by looking ahead the weight, unpromissing paths can be ignored in the 
early stage of the search. `atol` or `rtol` is used to control the accuracy 
for calculating the potential.

In addition, the algorithm assumes that WFST `T` is trimmed and the weight 
semiring is weakly left-divisible and zero-sum free such as tropical and 
log semirings.
"""
function push(T::WFST; atol::Real=0, rtol::Real= atol>0 ? 0 : 1e-3)
    return pushweight(T; atol, rtol)
end


