"""
    addstart(T::WFST{N,I,O,W}, s::N, w::V=one(W))
add start state `s` with weight `w` to `T`
"""
function addstart(T::WFST{N,I,O,W}, s::N, w::V=one(W)) where {N,I,O,W,V}
    T.starts[s] = w
    if !haskey(T.states, s)
        # states with no leaving arcs are still recorded with an empty
        # Arc Vector, this is important for some op like compose
        push!(T.states, s => Arc{N,I,O,W}[])
    end
    return nothing
end


"""
    addfinal(T::WFST{N,I,O,W}, f::N, w::V=one(V))
add final state `f` with weight `w` to `T`
"""
function addfinal(T::WFST{N,I,O,W}, f::N, w::V=one(W)) where {N,I,O,W,V}
    T.finals[f] = w
    if !haskey(T.states, f)
        # states with no leaving arcs are still recorded with an empty
        # Arc Vector, this is important for some op like compose
        push!(T.states, f => Arc{N,I,O,W}[])
    end
    return nothing
end



@inline rmstart(T::WFST{N}, k::N) where N = delete!(T.starts, k)
@inline rmfinal(T::WFST{N}, k::N) where N = delete!(T.finals, k)
@inline rmstate(T::WFST{N}, k::N) where N = delete!(T.states, k)

@doc """
    rmstart(T::WFST{N}, s::N)
remove start state `s` from `T`
""" rmstart

@doc """
    rmfinal(T::WFST{N}, f::N)
remove final state `f` from `T`
""" rmfinal

@doc """
    rmstate(T::WFST{N}, s::N)
remove state `s` and all its leaving arcs from `T`
""" rmstate



@inline function keyvecpush!(kvec::Dict{K, Vector{D}}, k::K, v::D) where {K,D}
    haskey(kvec, k) && 
    return push!(kvec[k], v );
    return push!(kvec,k=>[v]);
end


"""
    addeps(T::WFST{N,I,O,W}, src::N, dst::N, w::V=one(W))
add an epsilon arc from `src` state to `dst` state of `T` with weight `w`
"""
function addeps(T::WFST{N,I,O,W}, src::N, dst::N, w::V) where {N,I,O,W,V}
    keyvecpush!(T.states, src, Arc(dst, one(I), one(O), W(w)))
    !haskey(T.states, dst) && push!(T.states, dst => Arc{N,I,O,W}[])
    return nothing
end

function addeps(T::WFST{N,I,O,W}, src::N, dst::N, w::W=one(W)) where {N,I,O,W}
    keyvecpush!(T.states, src, Arc(dst, one(I), one(O), w))
    !haskey(T.states, dst) && push!(T.states, dst => Arc{N,I,O,W}[])
    return nothing
end


# add arc for an WFST from `Arc` format
"""
    addarc(T::WFST{N,I,O,W}, src::N, arc::Arc{N,I,O,V})
add leaving `arc` from `src` of `T`
"""
function addarc(T::WFST{N,I,O,W}, src::N, arc::Arc{N,I,O,V}) where {N,I,O,W,V}
    dst = next(arc)
    keyvecpush!(T.states, src, Arc(dst, ilabel(arc), olabel(arc), W(weight(arc))))
    !haskey(T.states, dst) && push!(T.states, dst => Arc{N,I,O,W}[])
    return nothing
end

function addarc(T::WFST{N,I,O,W}, src::N, arc::Arc{N,I,O,W}) where {N,I,O,W}
    dst = next(arc)
    keyvecpush!(T.states, src, arc)
    !haskey(T.states, dst) && push!(T.states, dst => Arc{N,I,O,W}[])
    return nothing
end


# add arc for an WFST from `src dst i o w` format
"""
    addarc(T::WFST{N,I,O,W}, src::N, dst::N, i::I, o::O, w::V=one(W))
add an arc from `src` to `dst` with weight `w`, input label `i` and output label `o` to `T`
"""
function addarc(T::WFST{N,I,O,W}, src::N, dst::N, i::I, o::O, w::V) where {N,I,O,W,V}
    keyvecpush!(T.states, src, Arc(dst, i, o, W(w)))
    !haskey(T.states, dst) && push!(T.states, dst => Arc{N,I,O,W}[])
    return nothing
end

function addarc(T::WFST{N,I,O,W}, src::N, dst::N, i::I, o::O, w::W=one(W)) where {N,I,O,W}
    keyvecpush!(T.states, src, Arc(dst, i, o, w))
    !haskey(T.states, dst) && push!(T.states, dst => Arc{N,I,O,W}[])
    return nothing
end


# add arc for an WFST from `src=>dst i o w` format
"""
    addarc(T::WFST{N,I,O,W}, src2dst::Pair{N,N}, i::I, o::O, w::V=one(W))
add an arc from `src` to `dst` with weight `w`, input label `i` and output label `o` to `T`, where (`src`, `dst`) is `src2dst`
"""
function addarc(T::WFST{N,I,O,W}, src2dst::Pair{N,N}, i::I, o::O, w::V=one(W)) where {N,I,O,W,V}
    src, dst = src2dst
    return addarc(T, src, dst, i, o, w)
end


# aid functions for parsing string into string
Base.parse(::Type{T}, str::String) where {T<:AbstractString} = T(str)
Base.parse(::Type{String}, str::String) = str
# add arc for an WFST from string of `src dst i:o/w` format
"""
    addarc(T::WFST, arcstr::AbstractString)
a convinient way to add an arc to `T` by string, mostly for quick prototyping
# Example
```julia
julia> begin
           T = WFST{Int,String,String,Float32}()
           addarc(T, "0 1 A:a/0.1") # 0 ── A:a/0.1 ──► 1
           addarc(T, "1 2 B:b")     # 1 ── B:b/1.0 ──► 2
           addarc(T, "2 3 C:/0.2")  # 2 ── C:ϵ/0.2 ──► 3
           addarc(T, "3 4 :d/0.3")  # 3 ── ϵ:d/0.3 ──► 4
           addarc(T, "4 5 :/0.4")   # 4 ── ϵ:ϵ/0.4 ──► 5
           addarc(T, "5 6 :")       # 5 ── ϵ:ϵ/1.0 ──► 6
           addarc(T, "6 7 E")       # 6 ── E:E/1.0 ──► 7
           addstart(T, 0)
           addfinal(T, 7)
       end
```
"""
function addarc(T::WFST{N,I,O,W}, arcstr::AbstractString) where {N,I,O,W}
    infos = string.( split(arcstr, " ", keepempty=false) )
    @assert length(infos)==3 begin
        "length of string-formed arc should contain ONLY 3 inputs"
    end
    p, n, arc = infos
    i, o, w = arc2iow(arc)
    src  = parse(N, p)
    dst  = parse(N, n)
    isym = parse(I, i)
    osym = parse(O, o)
    arcw = isequal(w, "") ? one(W) : parse(W, w)
    addarc(T, src, dst, isym, osym, arcw)
    return nothing
end


function arc2iow(arc::String)
    io_w = string.(split(arc, "/"))
    io   = first(io_w)
    i_o  = string.(split(io, ":"))
    i    = first(i_o)
    o = length( i_o)==2 ?  i_o[2] : i
    w = length(io_w)==2 ? io_w[2] : ""
    return i, o, w
end


"""
    addpath(T::WFST, path_linked_by_states_and_arcs::String)
Add a path to `T`. For example a path connected by states and arcs like \n
    [0]─a:A/0.7─►1─b:B/1─►2─c:c/0.9─►3─d:ϵ/0.7─►4─ϵ:E/1─►5─ϵ:ϵ/0.8─►6─ϵ:ϵ/1─►(7)
can be lazily added by command \n
    addpath(T, "0 a:A/0.7 1 b:B 2 c/0.9 3 d:/0.7 4 :E 5 :/0.8 6 : 7")
# Example
```julia
julia> begin
           T = WFST();
           addpath(T, "0 a:A/0.7 1 b:B 2 c/0.9 3 d:/0.7 4 :E 5 :/0.8 6 : 7");
           addstart(T,0); addfinal(T,7); T
       end
WFST{Float64}
 starts:
  0/1.0
 finals:
  7/1.0
 transitions:
  0→1 a:A/0.7
  1→2 b:B/1.0
  2→3 c:c/0.9
  3→4 d:ϵ/0.7
  4→5 ϵ:E/1.0
  5→6 ϵ:ϵ/0.8
  6→7 ϵ:ϵ/1.0
```
!!! note
    Arcs already added by `addpath` shall not be added again, i.e. if you are
    using `addpath` many times to construct a WFST, then paths shall not overlap 
    each other, otherwise the same arcs would appear many times.
"""
function addpath(T::WFST{N,I,O,W}, linked::String) where {N, I,O,W}
    path = string.(split(linked, " "))
    plen = length(path)
    @assert isodd(plen) && (plen > 1) begin
    "length of path info shall be like 3,5,7,9,11,..."
    end
    for k = 1 : 2 : plen-2
        p, arc, n = path[k:k+2]
        i, o, w = arc2iow(arc)
        src  = parse(N, p)
        dst  = parse(N, n)
        isym = parse(I, i)
        osym = parse(O, o)
        arcw = isequal(w, "") ? one(W) : parse(W, w)
        addarc(T, src, dst, isym, osym, arcw)
    end
end

