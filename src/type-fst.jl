
struct WFST{N,I,O,W}
    states :: Dict{N, Vector{Arc{N,I,O,W}}}
    starts :: Dict{N, W}
    finals :: Dict{N, W}
end


"""
    WFST{N,I,O,W}()
Make an `empty` WFST object, where \n
+ `N` is the type of states, e.g. `Int`
+ `I` is the type of input label on arc, e.g. `String`
+ `O` is the type of output label on arc, e.g. `String`
+ `W` is the type of weight on arc, e.g. `Float64`, `Semiring`
"""
function WFST{N,I,O,W}() where {N, I, O, W}
    states  = Dict{N, Vector{Arc{N,I,O,W}}}()
    startsw = Dict{N, W}()
    finalsw = Dict{N, W}()
    return WFST{N,I,O,W}(states, startsw, finalsw)
end


"""
    WFST(src::N, dst::N, i::I, o::O, w::W)
Make an `empty` WFST object via arg examples, where \n
+ `N` is the type of states, e.g. `Int`
+ `I` is the type of input label on arc, e.g. `String`
+ `O` is the type of output label on arc, e.g. `String`
+ `W` is the type of weight on arc, e.g. `Float64`, `Semiring`
!!! note
    all the input args are not used, they are just a type hint for convenience.
# Example
```julia
julia> WFST(0,1, "hello", "kitty", Float16(3.14))
WFST{Float16}
 starts:
 finals:
 transitions:
```
"""
function WFST(::N, ::N, ::I, ::O, ::W) where {N, I, O, W}
    states  = Dict{N, Vector{Arc{N,I,O,W}}}()
    startsw = Dict{N, W}()
    finalsw = Dict{N, W}()
    return WFST{N,I,O,W}(states, startsw, finalsw)
end


"""
    WFST()
Make an empty fst :: `WFST{Int,String,String,Float64}` object. This constructor is mainly used for fast prototyping.
"""
function WFST()
    states  = Dict{Int, Vector{Arc{Int,String,String,Float64}}}()
    startsw = Dict{Int, Float64}()
    finalsw = Dict{Int, Float64}()
    return WFST{Int,String,String,Float64}(states, startsw, finalsw)
end


"""
    WFST(states::Dict{N, Vector{Arc{N,I,O,W}}},
         starts::Dict{N, W},
         finals::Dict{N, W})

Make a fst object by passing all the field values. This constructor is mainly for `copy` and `deepcopy` use.

+ `N` is the type of states, e.g. `Int`
+ `I` is the type of input label on arc, e.g. `String`
+ `O` is the type of output label on arc, e.g. `String`
+ `W` is the type of weight on arc, e.g. `Float64`, `Semiring`
"""
function WFST(states::Dict{N, Vector{Arc{N,I,O,W}}},
              starts::Dict{N, W},
              finals::Dict{N, W}) where {N, I, O, W}
    return WFST{N,I,O,W}(states, starts, finals)
end




@inline cyan(  str::AbstractString) = "\e[36m" * str * "\e[0m"
@inline yellow(str::AbstractString) = "\e[33m" * str * "\e[0m"
@inline green( str::AbstractString) = "\e[32m" * str * "\e[0m"


function Base.show(io::IO, ::MIME"text/plain", f::WFST{N, I, O, W}) where {N, I, O, W}
    println(io, "WFST{$W}")
    println(io, yellow(" starts:"));i=1;
    for (state, weight) ∈ f.starts
        println(io, "  $state/$weight ")
        i += 1
        i > 10 && break
    end
    i > 10 && println("  ...")

    println(io, cyan(" finals:"));i=1;
    for (state, weight) ∈ f.finals
        println(io, "  $state/$weight ")
        i += 1
        i > 10 && break
    end
    i > 10 && println("  ...")

    println(io, green(" transitions:"));i=1;
    for (state, arcs) ∈ f.states
        j = 1
        for arc ∈ arcs
            println(io, "  $state", arc)
            j += 1
            j > 10 && break
        end
        i += 1
        i > 10 && break
    end
    i > 10 && println("  ...")
end


    Base.copy(x::WFST) = WFST(    copy(x.states),     copy(x.starts),     copy(x.finals))
Base.deepcopy(x::WFST) = WFST(deepcopy(x.states), deepcopy(x.starts), deepcopy(x.finals))


# judgements
@inline isstart(fsa::WFST{N}, s::N) where N = haskey(fsa.starts, s)
@inline isfinal(fsa::WFST{N}, f::N) where N = haskey(fsa.finals, f)

function haseps(T::WFST)
    for arcs ∈ values(T.states), arc ∈ arcs
        if isone(ilabel(arc)) || isone(olabel(arc))
            return true
        end
    end
    return false
end

function hasieps(T::WFST)
    for arcs ∈ values(T.states), arc ∈ arcs
        isone(ilabel(arc)) && return true
    end
    return false
end


function hasoeps(T::WFST)
    for arcs ∈ values(T.states), arc ∈ arcs
        isone(olabel(arc)) && return true
    end
    return false
end

function hasioeps(T::WFST)
    for arcs ∈ values(T.states), arc ∈ arcs
        if isone(ilabel(arc)) && isone(olabel(arc))
            return true
        end
    end
    return false
end


function iswfsa(fst::WFST) 
    for arcs ∈ values(T.states), arc ∈ arcs
        if ilabel(arc) ≠ olabel(arc)
            return false
        end
    end
    return true
end

# access
@inline starts(fsa::WFST) = fsa.starts
@inline finals(fsa::WFST) = fsa.finals

function ilabelsof(T::WFST{N, I, O, W}, keepeps::Bool=false) where {N, I, O, W}
    isyms = Set{I}()
    for arcs ∈ values(T.states)
        for arc ∈ arcs
            push!(isyms, ilabel(arc))
        end
    end
    !keepeps && setdiff!(isyms, [one(I)])
    return isyms
end


function olabelsof(T::WFST{N, I, O, W}, keepeps::Bool=false) where {N, I, O, W}
    osyms = Set{I}()
    for arcs ∈ values(T.states)
        for arc ∈ arcs
            push!(osyms, olabel(arc))
        end
    end
    !keepeps && setdiff!(osyms, [one(O)])
    return osyms
end


function iolabelsof(T::WFST{N, I, O, W}, keepeps::Bool=false) where {N, I, O, W}
    isyms = Set{I}()
    osyms = Set{I}()
    for arcs ∈ values(T.states)
        for arc ∈ arcs
            push!(isyms, ilabel(arc))
            push!(osyms, olabel(arc))
        end
    end
    !keepeps && setdiff!(isyms, [one(I)])
    !keepeps && setdiff!(osyms, [one(O)])
    return isyms, osyms
end


# numbers
@inline nstates(T::WFST) = length(T.states)
@inline nstarts(T::WFST) = length(T.starts)
@inline nfinals(T::WFST) = length(T.finals)
function narcs(T::WFST)
    n = 0
    for arcs ∈ values(T.states)
        n += length(arcs)
    end
    return n
end


"""
    length(T::WFST) -> nstates_plus_narcs :: Int
the sum of the number of states and arcs
"""
Base.length(T::WFST) = nstates(T) + narcs(T)


"""
    size(fst::WFST{N}) -> (nstates, narcs)
number of state and the number of arcs
"""
Base.size(T::WFST) = (nstates(T), narcs(T))


"""
    size(fst::WFST{N}, s::N)
number of arcs of state `s`
"""
Base.size(T::WFST{N}, s::N) where N = length(T.states[s])


function neps(T::WFST)
    c = 0
    for arcs ∈ values(T.states), arc ∈ arcs
        if isone(ilabel(arc)) || isone(olabel(arc))
            c += 1
        end
    end
    return c
end


function nieps(T::WFST)
    c = 0
    for arcs ∈ values(T.states), arc ∈ arcs
        if isone(ilabel(arc))
            c += 1
        end
    end
    return c
end

function noeps(T::WFST)
    c = 0
    for arcs ∈ values(T.states), arc ∈ arcs
        if isone(olabel(arc))
            c += 1
        end
    end
    return c
end


function nioeps(T::WFST)
    c = 0
    for arcs ∈ values(T.states), arc ∈ arcs
        if isone(ilabel(arc)) && isone(olabel(arc))
            c += 1
        end
    end
    return c
end

