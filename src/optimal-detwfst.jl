"""
    detwfst(T::WFST{N, I, O, W}) -> D::WFST{Vector{Tuple{N,O,W}}, I, O, W}

Determinize `T`, producing a new deterministic **WFST** `D`. A deterministic WFST has **ONLY** one path 
for a given input sequence. In a nondeterministic WFST, a state with multiple leaving arcs having 
the same input label introduces ambiguity. The determinization process eliminates this ambiguity 
by ensuring that each input symbol corresponds to a unique leaving arc at any state. After `detwfst`, 
calling `reorder` on `D` is recommoded to re-map the state to integer indices. It has an alias `tdet`.

!!! note
    + If `T` contains ϵ:ϵ-labeled transitions, ensure they are removed before calling `detwfst`.
    + `ϵ` labels are treated as regular symbols.
    + The transducer must be functional.
    + The weights must be weakly left divisible (valid for Tropical and Log Semirings)
    + The weights must be zero-sum-free, i.e. if ```x ⊕ y = 0````, then ```x = y = 0```.
    + The input WFST `T` is not p-subsequential.

# Design
    ───────────────────────────────────────────────────────────────────────
    pzvs = {(p,z,v)}, a set of weighted src state with output z
    qzvs = {(q,z,v)}, a set of weighted dst state with output z
    where    p denotes source state,
             q denotes destination state, 
             z denotes leftover output label and
             v denotes the residual weight of p or q.
    ───────────────────────────────────────────────────────────────────────
    contains two main loops:
        ● 1-st level: iteration of the weighted set of states (p,z,v) ∈ Q[p′]
        ● 2-nd level: e ∈ E(p)
    ───────────────────────────────────────────────────────────────────────
"""
function detwfst(T::WFST{N,I,O,W}) where {N,I,O,W}
    Ṅ = Vector{Tuple{N,O,W}}
    T = onestart(T)
    D = WFST{Ṅ, I, O, W}()
    Q = Vector{Ṅ}()  # queue for doing visiting, Vector struct is fine
    S = Set{Ṅ}()     # set of weighted state Set for existence checking
    ε = one(O)
    for (p, v) ∈ T.starts
        pzv = [(p, ε, v)]
        push!(Q, pzv)
        addstart(D, pzv)
    end

    E = T.states  # state => edges
    F = T.finals  # state => weight
    OW = Tuple{O,W}
    while !isempty(Q)
        pzvs = pop!(Q)
        xyw  = Dict{I, OW}()
        xqyw = Dict{I, Dict{N,OW}}()
        for (p, z, v) ∈ pzvs, e ∈ E[p]
            x = ilabel(e)
            y = olabel(e)
            w = weight(e)
            q  =  next(e)
            ỹw̌ = tuple(z * y, v * w)
            # ⊕ x:y/w info
            if haskey(xyw, x)
                xyw[x] = xyw[x] .+ ỹw̌
            else
                xyw[x] = ỹw̌
            end
            # ⊕ x:y/w ─► q info
            if haskey(xqyw, x)
                if haskey(xqyw[x], q)
                    xqyw[x][q] = xqyw[x][q] .+ ỹw̌
                else
                    xqyw[x][q] = ỹw̌
                end
            else
                xqyw[x] = Dict{N,OW}(q => ỹw̌)
            end
        end

        for (x, (y,w)) ∈ xyw
            # qzvs holds residual weight and leftover output (would be used as prefix)
            qzvs = Vector{Tuple{N,O,W}}()
            for (q, (ỹ,w̌)) ∈ xqyw[x]
                #= residual weight w⁻¹ * w̌, where w is the `prefix` weight but w⁻¹ is not
                pre-computed, because its type may be an inexpressible type, e.g String =#
                push!(qzvs, (q, ỹ/y, w̌/w))
            end
            addarc(D, pzvs, qzvs, x, y, w)

            if qzvs ∉ S
                push!(S, qzvs)
                push!(Q, qzvs)
            end

            hasfinal = false
            finalval = zero(W)
            for (q, z, v) ∈ qzvs
                if isfinal(T, q)
                    hasfinal = true
                    finalval += v * F[q]
                end
            end
            hasfinal && addfinal(D, qzvs, finalval)
        end
    end
    return D
end



@inline function Base.:+(x::AbstractString, y::AbstractString)
    !isequal(x, y) && return ""
    return x
end


# @inline function Base.:/(x::AbstractString, y::AbstractString)
#     !isequal(x, y) && return x
#     return ""
# end


# function Base.:+(x::AbstractString, y::AbstractString)
#     offset = 0  # offset from 1st index
#     n = min(length(x), length(y))
#     for i = 1:n
#         !isequal(x[i], y[i]) && break
#         offset = i
#     end
#     return x[1 : offset]
# end

function Base.:/(x::AbstractString, y::AbstractString)
    lx = length(x)
    ly = length(y)
    ly > lx && error("length of \"$x\" shall be greater than length of \"$y\"")

    offset = 0  # offset from 1st index
    for i = 1:ly
        !isequal(x[i], y[i]) && break
        offset = i
    end
    return x[offset+1:lx]
end



"""
    tdet(T::WFST)
`tdet` is an alias of `detwfst`
"""
tdet(T::WFST) = detwfst(T)


function det(T::WFST{N,I,O,W}) where {N,I,O,W}
    Ṅ = Vector{Tuple{N,O,W}}
    T = onestart(T)
    D = WFST{Ṅ, I, O, W}()
    Q = Vector{Ṅ}()  # queue for doing visiting, Vector struct is fine
    S = Set{Ṅ}()     # set of weighted state Set for existence checking
    ε = one(O)
    for (p, v) ∈ T.starts
        pzv = [(p, ε, v)]
        push!(Q, pzv)
        addstart(D, pzv)
    end
    M = maxid(T)
    E = T.states  # state => edges
    F = T.finals  # state => weight
    NO = Tuple{N,O}
    OW = Tuple{O,W}
    while !isempty(Q)
        pzvs = pop!(Q)
        xyw  = Dict{I, OW}()
        xqyw = Dict{I, Dict{NO,W}}()
        for (p, z, v) ∈ pzvs, e ∈ E[p]
            q =   next(e)
            x = ilabel(e)
            y = olabel(e)
            w = weight(e)
            w̌ = v * w

            # ⊕ x:y/w info
            ỹw̌ = isone(z) ? tuple(y, w̌) : tuple(z, w̌)
            if haskey(xyw, x)
                xyw[x] = xyw[x] .+ ỹw̌
            else
                xyw[x] = ỹw̌
            end
            # ⊕ x:y/w ─► q info
            qy = tuple(q, z * y)
            if haskey(xqyw, x)
                if haskey(xqyw[x], qy)
                    xqyw[x][qy] = xqyw[x][qy] + w̌
                else
                    xqyw[x][qy] = w̌
                end
            else
                xqyw[x] = Dict{NO,W}(qy => w̌)
            end
        end

        for (x, (y,w)) ∈ xyw
            # qzvs is dst states' container with residual weight and leftover output (would be used as prefix)
            qzvs = Vector{Tuple{N,O,W}}()
            for ((q,ỹ),w̌) ∈ xqyw[x]
                #= residual weight w⁻¹ * w̌, where w is the `prefix` weight but w⁻¹ is not
                pre-computed, because its type may be an inexpressible type, e.g String =#
                v = !isfinal(T, q) ? (w̌/w) : (w̌/w * F[q])
                push!(qzvs, (q, ỹ/y, v))
            end
            addarc(D, pzvs, qzvs, x, y, w)

            if qzvs ∉ S
                push!(S, qzvs)
                push!(Q, qzvs)
            end

            for (q, z, v) ∈ qzvs
                if isfinal(T, q)
                    addfinal(D, qzvs)
                    break
                end
            end
        end
    end

    oldfinal = []
    newfinal = []
    for pzvs ∈ keys(D.finals), (p, z, v) ∈ pzvs
        !isfinal(T, p) && continue
        isone(z) && continue
        M = M + 1
        qzvs = [(M,one(O),one(W))]
        addarc(D, pzvs, qzvs, one(I), z, v)
        push!(newfinal, qzvs)
        push!(oldfinal, pzvs)
    end

    for pzvs ∈ oldfinal rmfinal(D, pzvs) end
    for pzvs ∈ newfinal addfinal(D, pzvs) end

    return D
end

