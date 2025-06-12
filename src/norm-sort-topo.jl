"""
    indegree(T::WFST{N})
Return the indegree of all states in `T`
"""
function indegree(T::WFST{N}) where N
    E  = T.states
    qi = Dict{N,Int}()  # Dict: state => indegrees
    qs = Vector{N}()    # state container for iterating

    for start ∈ keys(T.starts)
        push!(qs, start)
    end

    while !isempty(qs)
        p = pop!(qs)    # or popfirst!
        for e ∈ E[p]
            q = next(e)
            if !haskey(qi, q)
                push!(qs, q)    # q was not in the queue, ensure enqueue ONLY once
                qi[q] = 1       # maybe q has only one incoming ϵ-arc
            else
                qi[q] += 1      # here q has multiple incoming ϵ-arcs
            end
        end
    end
    return qi
end


"""
    toposort(T::WFST{N}) -> sorted::Vector{N}
Return states of `T` in topological order. `T` must be acyclic.
"""
function toposort(T::WFST{N}) where N
    E  = T.states
    qi = indegree(T)        # Dict: state => indegrees
    qs = Vector{N}()        # state container for iterating
    rank = Vector{N}()      # topological ordered states

    for start ∈ keys(T.starts)
        push!(qs, start)
    end

    while !isempty(qs)
        p = pop!(qs)   # chose popfirst! or pop!
        push!(rank, p) # p is pushed right after poping
        for e ∈ E[p]
            q = next(e)
            qi[q] -= 1
            if iszero(qi[q])
                push!(qs, q)
            end
        end
    end
    return rank
end


"""
    levelsort(T::WFST{N}) -> sorted::Vector{Vector{N}}
Return level sorted states of `T` in topological order. `T` must be acyclic.
"""
function levelsort(T::WFST{N}) where N
    E  = T.states
    qi = indegree(T)            # Dict: state => indegrees
    qs = Vector{Vector{N}}()    # state container for iterating
    rank = Vector{Vector{N}}()  # state Vector container for iterating

    inits = Vector{N}()
    for s ∈ keys(T.starts)
        push!(inits, s)
    end

    push!(qs, inits)
    while !isempty(qs)
        ps = pop!(qs)
        push!(rank, ps)
        level = Vector{N}()
        for p ∈ ps
            for e ∈ E[p]
                q = next(e)
                qi[q] -= 1
                if iszero(qi[q])
                    push!(level, q)
                end
            end
        end
        !isempty(level) && push!(qs, level)
    end
    return rank
end

