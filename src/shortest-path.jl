"""
A linked list tracer for state and output label
"""
mutable struct Tracer{N,O}
    state :: N  # state
    label :: O  # olabel
    previous :: Union{Tracer, Nothing}
    function Tracer(p::N, o::O) where {N,O}
        new{N,O}(p, o, nothing)
    end
end

@inline stateof(t::Tracer)  = t.state
@inline labelof(t::Tracer)  = t.label
@inline previous(t::Tracer) = t.previous

function addprev(t::Tracer{N,O}, prev::Tracer{N,O}) where {N,O}
    t.previous = prev
end

function Base.show(io::IO, t::Tracer)
    p = stateof(t)
    o = labelof(t)
    print(io, "Tracer($p, $o)")
end

function printpath(t::Tracer)
    olabels = Vector{String}()
    push!(olabels, labelof(t))
    p = previous(t)
    while !isnothing(p)
        push!(olabels, labelof(p))
        p = previous(p)
    end
    str = ""
    for s ∈ reverse(olabels)
        str *= s
        str *= " "
    end
    println(str)
end


"""
priority queue for k-best-strings problem.
"""
mutable struct PathQueue{N,O,W}
    cmper :: Function
    dists :: Dict{N,W}
    queue :: Vector{Tuple{Tracer{N,O},W}}
    function PathQueue(T::WFST{N,I,O,W}, cmpf::Function = <,
                                         atol::Real = 0, 
                                         rtol::Real = atol>0 ? 0 : 1e-3) where {N,I,O,W}
        δists = potential(T; atol, rtol)
        queue = Vector{Tuple{Tracer{N,O},W}}()
        for (start, w) ∈ T.starts
            tracer = Tracer(start, one(O))
            push!(queue, (tracer, w))
        end
        new{N,O,W}(cmpf, δists, queue)
    end
end

@inline Base.length(P::PathQueue)  = length(P.queue)
@inline Base.isempty(P::PathQueue) = isempty(P.queue)

function Base.push!(P::PathQueue{N,O,W}, no_w::Tuple{Tracer{N,O},W}) where {N,O,W}
    return push!(P.queue, no_w)
end

# TODO: a faster priority method for large k is needed
function Base.pop!(P::PathQueue{N,O,W}) where {N,O,W}
    f = P.cmper
    V = P.dists
    Q = P.queue
    INDX = zero(Int)
    LAST = length(Q)     #   worst -∞ : ∞  best
    BEST = isequal(f, >) ? typemin(W) : typemax(W)
    for (i, (tracer,c)) ∈ enumerate(Q)
        p = stateof(tracer)
        score = c * V[p]
        if f(score, BEST) # f(x,y)==true ⇒ x is better than y
            BEST = score
            INDX = i
        end
    end
    Q[LAST], Q[INDX] = 
    Q[INDX], Q[LAST];    
    return pop!(Q)
end


"""
    paths(T::WFST{N,I,O,W},
          k::Int;
          cmpf::Function = <,
          atol::Real = 0,
          rtol::Real = atol>0 ? 0 : 1e-3) -> wtracers::Vector{ Tuple{Tracer{N,O}, W} }

Return the weighted `k` best strings. 
If `cmpf` is `>`, then  bigger is better (e.g. prob Semiring), 
if `cmpf` is `<`, then smaller is better (e.g. min-plus Semiring)
# Example
```julia
julia> am  = amfst(3, 6);
julia> wtr = paths(am, 3, cmpf= >);
julia> for (tr,w) ∈ wtr
            print("weight=\$w, path:")
            printpath(tr)
       end
weight=0.023901694, path: s₂ s₃ s₁ s₂ s₃ s₁
weight=0.020913983, path: s₂ s₃ s₃ s₂ s₃ s₁
weight=0.019121354, path: s₂ s₃ s₁ s₂ s₃ s₂
```
"""
function paths(T::WFST{N,I,O,W}, k::Int; cmpf::Function = <,
                                         atol::Real = 0, 
                                         rtol::Real = atol>0 ? 0 : 1e-3) where {N,I,O,W}
    #= r[q] is counts since last time q was extracted from S
    Q is a priority queue that pops best (po,c) out every time =#
    r = Dict{N,Int}(q=>0 for q ∈ keys(T.states))
    Q = PathQueue(T, cmpf, atol, rtol)
    wtracers = Vector{ Tuple{Tracer{N,O}, W} }()
    ftouched = 0  # count of touched final states
    E = T.states
    F = T.finals    
    while !isempty(Q)
        tracer, c = pop!(Q)
        p = stateof(tracer)
        #= every state shall only be went throught at most k times,
        i.e. each path of the k-best-paths goes throught it at most once. so, 
        [1] terminate if any final state is extended k times and
        [2] don't extend if any state is extended more than k times =#
        r[p] = r[p] + 1
        r[p] > k && continue
        if isfinal(T, p)
            ftouched += 1
            push!(wtracers, (tracer, c * F[p]))
            isequal(ftouched, k) && break;
        end
        # extend if r[p] ≤ k
        for e ∈ E[p]
            n = next(e)
            c̄ = c * weight(e)
            head = Tracer(n, olabel(e))
            addprev(head, tracer)
            push!(Q, (head, c̄))
        end
    end
    return wtracers
end




# demo usage
function sumnorm(x::Array)
    return x .* inv(sum(x))
end

function amfst(N::Int, T::Int)
    A = WFST{Int,String,String,Float32}()
    addstart(A, 0)
    addfinal(A, T)
    for t = 1:T
        i = "x" * subdigits(t)
        v = sumnorm(rand(Float32, N))
        for n = 1:N
            o = "s" * subdigits(n)
            w = floor(v[n], digits=2)
            addarc(A, t-1, t, i, o, w)
        end
    end
    return A
end

function hmmfst(N::Int)
    A = WFST{Int,String,String,Float32}()
    for n = 1:N
        i = "s" * subdigits(n)
        v = sumnorm(rand(Float32, N))
        for m = 1:N
            o = "s" * subdigits(m)
            w = floor(v[m], digits=2)
            addarc(A, n, m, i, o, w)
        end
        addstart(A, n)
        addfinal(A, n)
    end
    return A
end
