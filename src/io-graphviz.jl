SUB = Dict{Char, String}(
    '0'=>"₀",
    '1'=>"₁",
    '2'=>"₂",
    '3'=>"₃",
    '4'=>"₄",
    '5'=>"₅",
    '6'=>"₆",
    '7'=>"₇",
    '8'=>"₈",
    '9'=>"₉"
);


"""
    subdigits(n::Int)::String
Convert the decimal number to lower right corner digits in the form of string

# Example
    julia> "X" * subdigits(321)
    "X₃₂₁"
"""
function subdigits(n::Integer)
    s = ""
    for c ∈ string(n)
        s *= SUB[c]
    end
    return s
end



function tupletxt(etc::T) where T <:Tuple
    t = "("
    n = length(etc)
    for e ∈ etc[1:n-1]
        t *= "$e,"
    end;t *= "$(etc[n]))"
    return t
end

@inline ᵖ(state::Int) = "$state"

function ᵖ(state::Vector{T}) where T <:Tuple
    s = "["
    n = length(state)
    for e ∈ state[1:n-1]
        s *= "$(tupletxt(e)),"
    end;s *= "$(tupletxt(state[n]))]"
    return s
end


"""
    dot(wfst::WFST; rankdir::String="LR",
                    head::String="digraph",
                    name::String="WFST",
                    keepone ::Bool=true,
                    keepsame::Bool=true,
                    offset::Int=0,
                    digits::Int=2,
                    edgefont::String="Arial",
                    nodefont::String="Arial")

convert a WFST into dot format string

# Arguments
`rankdir` :
+ `"LR"` is left to right layout → \n
+ `"RL"` is right to left layout ← \n
+ `"TB"` is top to down layout ↓ \n
+ `"BT"` is down to top layout ↑ \n
`head` : has two options, "digraph" or "subgraph"\n
`name` : label for this graph\n
`keepsame` : if `false`, then only show one symbol on arc when i o labels are the same.\n
`offset` : ignore this, for dev only \n
`edgefont` : font name of edge, like Arial/Georgia/微软雅黑/宋体/楷体/仿宋/Aldhabi etc \n
`nodefont` : font name of node, see `edgefont`, just use the one you like \n

# Online visualization with Graphviz
there are some websides providing free online preview service:
+ `https://edotor.net/` (highly recommended)

+ `https://dreampuf.github.io/GraphvizOnline/` (normal)
+ `https://obren.io/tools/graphviz/` (normal)

+ `https://www.devtoolsdaily.com/graphviz/`  (enhanced coloring with AI)
+ `http://magjac.com/graphviz-visual-editor/` (interactive and editable)

+ `https://sketchviz.com/new` (sketch style)
"""
function dot(T::WFST{N,IL,OL,W}; rankdir::String="LR",
                                 head::String="digraph",
                                 name::String="",
                                 keepone ::Bool=true,
                                 keepsame::Bool=true,
                                 offset::Int=0,
                                 digits::Int=2,
                                 startcolor::String="gold",
                                 finalcolor::String="pink",
                                 edgefont::String="Arial",
                                 nodefont::String="Arial") where {N,IL,OL,W}
    I = T.starts
    F = T.finals
    s = ""

    bias = if iszero(offset)
        ""
    else
        subdigits(offset)
    end

    if head == "subgraph"
        head *= " cluster_$offset"
        s *= "$head \n{\n\tlabel = \"$name\";\n"
    else
        s *= "$head \n{\n\trankdir = $rankdir;\n\tbgcolor = \"transparent\";\n\tlabel = \"$name\";\n"
    end

    s *= "\tedge [arrowsize=0.7, fontname=\"$edgefont\"]\n"
    s *= "\tnode [color=black, fontname=\"$nodefont\"]\n"

    # dealing with arcs leaving src
    for (src, arcs) ∈ T.states
        # handle starts and finals
        bestart = isstart(T, src)
        befinal = isfinal(T, src)
        psrc = ᵖ(src)
        ID = psrc * "$bias"
        s *= "\t\"$(ID)\" "
        if bestart && befinal
            wi = round(I[src]; digits)
            wf = round(F[src]; digits)
            if !isone(wi) || !isone(wf)
                s *= "[style=filled,fillcolor=orangered,shape=doublecircle, penwidth=2.5, label=\"$psrc/$wi,$wf\"];\n"
            else
                # wi & wf are all one, then don't show weights
                s *= "[style=filled,fillcolor=orangered,shape=doublecircle, penwidth=2.5, label=\"$psrc\"];\n"
            end
        elseif bestart
            w = round(I[src]; digits)
            v = ifelse(isone(w) && !keepone, "" , "/$w")
            s *= "[style=filled,fillcolor=$startcolor,shape=circle, penwidth=2.5, label=\"$psrc$v\"];\n"
        elseif befinal
            w = round(F[src]; digits)
            v = ifelse(isone(w) && !keepone, "" , "/$w")
            s *= "[style=filled,fillcolor=$finalcolor,shape=doublecircle, label=\"$psrc$v\"];\n"
        else
            s *= "[shape=circle, label=\"$psrc\"];\n"
        end

        for arc ∈ arcs
            dst = ᵖ(next(arc)) * "$bias"
            i = ilabel(arc)
            o = olabel(arc)
            w = weight(arc)
            w = round(w; digits)
            is = ifelse(isone(i), "ε", i)
            os = ifelse(isone(o), "ε", o)
            ws = ifelse(isone(w) && !keepone, "" , "/$w")

            e = ifelse(is==os && !keepsame, "$is$ws", "$is:$os$ws")
            s = "$s\t\t\"$ID\" -> \"$dst\" [label=\"$e\"];\n"
        end
    end

    return s * "}\n"
end


function savedot(dotname::String,
                 T::WFST{N,I,O,W};
                 rankdir::String="LR",
                 head::String="digraph",
                 name::String="",
                 keepone ::Bool=true,
                 keepsame::Bool=true,
                 offset::Int=0,
                 digits::Int=2,
                 edgefont::String="Arial",
                 nodefont::String="Arial") where {N,I,O,W}
    body = dot(T; rankdir, head, name, keepone, keepsame, offset, digits, edgefont, nodefont)
    open(dotname, "w") do io
        write(io, body)
    end
    return nothing
end


"""
    dots(wfsts::Vector, names::Vector{String}; rankdir::String="LR",
                                               digits::Int=2,
                                               keepsame::Bool=true,
                                               keepone ::Bool=true,
                                               edgefont::String="Arial",
                                               nodefont::String="Arial")

convert several WFSTs into dot format string

# Arguments
`rankdir` :
+ "LR" is left to right layout → \n
+ "RL" is right to left layout ← \n
+ "TB" is top to down layout ↓ \n
+ "BT" is down to top layout ↑ \n
`names` : labels for each sub graph\n
`keepsame` : if `false`, then only show one symbol on arc when i o labels are the same.\n

!!! note
    We can NOT control each subgraph's direction independently. `direction` only controls global orientation.
"""
function dots(wfsts::Vector,
              names::Vector{String};
              rankdir::String="LR",
              keepsame::Bool=true,
              keepone::Bool=true,
              edgefont::String="Arial",
              nodefont::String="Arial",
              digits::Int=2)
    s = "digraph WFSTs \n{\nrankdir = $rankdir;\nbgcolor = \"transparent\";\n\n"
    n = 0
    for (i, wfst) ∈ enumerate(wfsts)
        s *= dot(wfst, rankdir=rankdir,
                       head="subgraph",
                       name=names[i],
                       keepsame=keepsame,
                       keepone=keepone,
                       offset=n,
                       edgefont=edgefont,
                       nodefont=nodefont,
                       digits=digits)
        n += nstates(wfst)
    end
    s = s * "}\n"
    return s
end


function savedots(dotname::String, names::Vector{String}, T::WFST{N,I,O,W}; kws...) where {N,I,O,W}
    body = dots(T, names; kws...)
    open(dotname, "w") do io
        write(io, body)
    end
    return nothing
end


"""
    save(T::WFST, filename::String; direction::String="LR")

Draw a WFST into `filename`, format info shall be included in filename. Graphviz shall be intalled.
# Arguments
`filename`: file's name, like `"./home/data/xxx.svg"`
"""
function savewfst(filename::String, T::WFST; rankdir::String="LR")
    dotstring = dot(T; rankdir)
    fmt = last(split(filename, "."))
    isone(fmt) && error("the name must have a format suffix, e.g. pdf/svg")
    run(`echo $dotstring` |> `dot -T $fmt -o $filename`)
end




function leveldot(T::WFST{N,IL,OL,W};
                  rankdir::String="LR",
                  head::String="digraph",
                  name::String="",
                  keepone ::Bool=true,
                  keepsame::Bool=true,
                  offset::Int=0,
                  digits::Int=2,
                  startcolor::String="gold",
                  finalcolor::String="pink",
                  edgefont::String="Arial",
                  nodefont::String="Arial") where {N,IL,OL,W}
    I = T.starts
    F = T.finals
    E = T.states

    s = ""

    bias = subdigits(offset)
    if head == "subgraph"
        head *= " cluster_$offset"
        s *= "$head \n{\n\tlabel = \"$name\";\n"
    else
        s *= "$head \n{\n\trankdir = $rankdir;\n\tbgcolor = \"transparent\";\n\tlabel = \"$name\";\n"
    end

    s *= "\tedge [arrowsize=0.7, fontname=\"$edgefont\"]\n"
    s *= "\tnode [color=black, fontname=\"$nodefont\"]\n"

    # dealing with arcs leaving src
    for srcvec ∈ levelsort(T)
        ranked = "\t{rank=\"same\""
        for src ∈ srcvec
            # handle starts and finals
            bestart = isstart(T, src)
            befinal = isfinal(T, src)

            ID = "$src$bias"
            s *= "\t\"$(ID)\" "
            if bestart && befinal
                wi = round(I[src]; digits)
                wf = round(F[src]; digits)
                vi = ifelse(isone(wi) && !keepone, "" , "$wi")
                vf = ifelse(isone(wf) && !keepone, "" , "$wf")
                s *= "[style=filled,fillcolor=orangered,shape=doublecircle, penwidth=2.5, label=\"$src/$vi,$vf\"];\n"
            elseif bestart
                w = round(I[src]; digits)
                v = ifelse(isone(w) && !keepone, "" , "/$w")
                s *= "[style=filled,fillcolor=$startcolor,shape=circle, penwidth=2.5, label=\"$src$v\"];\n"
            elseif befinal
                w = round(F[src]; digits)
                v = ifelse(isone(w) && !keepone, "" , "/$w")
                s *= "[style=filled,fillcolor=$finalcolor,shape=doublecircle, label=\"$src$v\"];\n"
            else
                s *= "[shape=circle, label=\"$src\"];\n"
            end

            for arc ∈ E[src]
                dst = string(next(arc)) * "$bias"
                i = ilabel(arc)
                o = olabel(arc)
                w = weight(arc)
                w = round(w; digits)
                is = ifelse(isone(i), "ε", i)
                os = ifelse(isone(o), "ε", o)
                ws = ifelse(isone(w) && !keepone, "" , "/$w")

                e = ifelse(is==os && !keepsame, "$is$ws", "$is:$os$ws")
                s = "$s\t\t\"$ID\" -> \"$dst\" [label=\"$e\"];\n"
            end
            ranked *= "; \"$ID\""
        end
        ranked *= "}\n"
        s *= ranked
    end

    return s * "}\n"
end


function saveleveldot(dotname::String,
                 T::WFST{N,I,O,W};
                 rankdir::String="LR",
                 head::String="digraph",
                 name::String="",
                 keepone ::Bool=true,
                 keepsame::Bool=true,
                 offset::Int=0,
                 digits::Int=2,
                 edgefont::String="Arial",
                 nodefont::String="Arial") where {N,I,O,W}
    body = leveldot(T; rankdir, head, name, keepone, keepsame, offset, digits, edgefont, nodefont)
    open(dotname, "w") do io
        write(io, body)
    end
    return nothing
end


"""
    all states are placed in one line.
"""
function linedot(T::WFST{N,IL,OL,W}; rankdir::String="LR",
                                 head::String="digraph",
                                 name::String="",
                                 keepone ::Bool=true,
                                 keepsame::Bool=true,
                                 offset::Int=0,
                                 digits::Int=2,
                                 edgefont::String="Arial",
                                 nodefont::String="Arial") where {N,IL,OL,W}
    I = T.starts
    F = T.finals
    s = ""

    bias = if iszero(offset)
        ""
    else
        subdigits(offset)
    end

    ranked = "\t{rank=\"same\""

    if head == "subgraph"
        head *= " cluster_$offset"
        s *= "$head \n{\n\tlabel = \"$name\";\n"
    else
        s *= "$head \n{\n\trankdir = $rankdir;\n\tbgcolor = \"transparent\";\n\tlabel = \"$name\";\n"
    end

    s *= "\tedge [arrowsize=0.7, fontname=\"$edgefont\"]\n"
    s *= "\tnode [color=black, fontname=\"$nodefont\"]\n"

    # dealing with arcs leaving src
    for (src, arcs) ∈ T.states
        # handle starts and finals
        bestart = isstart(T, src)
        befinal = isfinal(T, src)

        ID = "$src$bias"
        s *= "\t\"$(ID)\" "
        ranked *= "; \"$ID\""

        if bestart && befinal
            wi = round(I[src]; digits)
            wf = round(F[src]; digits)
            if !isone(wi) || !isone(wf)
                s *= "[shape=doublecircle, penwidth=2.5, label=\"$src/$wi,$wf\"];\n"
            else
                # wi & wf are all one, then don't show weights
                s *= "[shape=doublecircle, penwidth=2.5, label=\"$src\"];\n"
            end
        elseif bestart
            w = round(I[src]; digits)
            v = ifelse(isone(w) && !keepone, "" , "/$w")
            s *= "[shape=circle, penwidth=2.5, label=\"$src$v\"];\n"
        elseif befinal
            w = round(F[src]; digits)
            v = ifelse(isone(w) && !keepone, "" , "/$w")
            s *= "[shape=doublecircle, label=\"$src$v\"];\n"
        else
            s *= "[shape=circle, label=\"$src\"];\n"
        end

        for arc ∈ arcs
            dst = string(next(arc)) * "$bias"
            i = ilabel(arc)
            o = olabel(arc)
            w = weight(arc)
            w = round(w; digits)
            is = ifelse(isone(i), "ε", i)
            os = ifelse(isone(o), "ε", o)
            ws = ifelse(isone(w) && !keepone, "" , "/$w")

            e = ifelse(is==os && !keepsame, "$is$ws", "$is:$os$ws")
            s = "$s\t\t\"$ID\" -> \"$dst\" [label=\"$e\"];\n"
        end
    end
    ranked *= "}\n"
    s *= ranked
    return s * "}\n"
end


function depthdot(T::WFST{N,IL,OL,W},
                  direction::String;
                  rankdir::String="LR",
                  head::String="digraph",
                  name::String="",
                  keepone ::Bool=true,
                  keepsame::Bool=true,
                  offset::Int=0,
                  digits::Int=2,
                  startcolor::String="gold",
                  finalcolor::String="pink",
                  edgefont::String="Arial",
                  nodefont::String="Arial") where {N,IL,OL,W}
    I = T.starts
    F = T.finals
    E = T.states

    s = ""

    bias = subdigits(offset)
    if head == "subgraph"
        head *= " cluster_$offset"
        s *= "$head \n{\n\tlabel = \"$name\";\n"
    else
        s *= "$head \n{\n\trankdir = $rankdir;\n\tbgcolor = \"transparent\";\n\tlabel = \"$name\";\n"
    end

    s *= "\tedge [arrowsize=0.7, fontname=\"$edgefont\"]\n"
    s *= "\tnode [color=black, fontname=\"$nodefont\"]\n"

    depth2states, maxdepth = depthvecs(T, direction)
    depthiter = if isequal(direction, "fwd")
        (0 : +1 : maxdepth)
    else
        (maxdepth : -1 : 0)
    end

    # dealing with arcs leaving src
    for depth ∈ depthiter
        srcvec =  depth2states[depth]
        ranked = "\t{rank=\"same\""
        for src ∈ srcvec
            # handle starts and finals
            bestart = isstart(T, src)
            befinal = isfinal(T, src)

            ID = "$src$bias"
            s *= "\t\"$(ID)\" "
            if bestart && befinal
                wi = round(I[src]; digits)
                wf = round(F[src]; digits)
                vi = ifelse(isone(wi) && !keepone, "" , "$wi")
                vf = ifelse(isone(wf) && !keepone, "" , "$wf")
                s *= "[style=filled,fillcolor=orangered,shape=doublecircle, penwidth=2.5, label=\"$src/$vi,$vf\"];\n"
            elseif bestart
                w = round(I[src]; digits)
                v = ifelse(isone(w) && !keepone, "" , "/$w")
                s *= "[style=filled,fillcolor=$startcolor,shape=circle, penwidth=2.5, label=\"$src$v\"];\n"
            elseif befinal
                w = round(F[src]; digits)
                v = ifelse(isone(w) && !keepone, "" , "/$w")
                s *= "[style=filled,fillcolor=$finalcolor,shape=doublecircle, label=\"$src$v\"];\n"
            else
                s *= "[shape=circle, label=\"$src\"];\n"
            end

            for arc ∈ E[src]
                dst = string(next(arc)) * "$bias"
                i = ilabel(arc)
                o = olabel(arc)
                w = weight(arc)
                w = round(w; digits)
                is = ifelse(isone(i), "ε", i)
                os = ifelse(isone(o), "ε", o)
                ws = ifelse(isone(w) && !keepone, "" , "/$w")

                e = ifelse(is==os && !keepsame, "$is$ws", "$is:$os$ws")
                s = "$s\t\t\"$ID\" -> \"$dst\" [label=\"$e\"];\n"
            end
            ranked *= "; \"$ID\""
        end
        ranked *= "}\n"
        s *= ranked
    end

    return s * "}\n"
end
