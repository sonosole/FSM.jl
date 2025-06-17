"""
    mermaid(T::WFST{N,IL,OL,W}; rankdir::String="LR",
                                keepone ::Bool=true,
                                keepsame::Bool=true,
                                digits::Int=2) where {N,IL,OL,W}

Convert a WFST into mermaid flow graph string code, which could 
be embedded in Markdown file.

# Keyword Arguments
`rankdir` :
+ `"LR"` is left to right layout → \n
+ `"RL"` is right to left layout ← \n
+ `"TB"` is top to down layout ↓ \n
+ `"BT"` is down to top layout ↑ \n
`keepsame` : if `false`, then only show one symbol on arc when i o labels are the same.\n
`digits` : number of digits kept for the weights afer decimal point, e.g. if digits=2, then 3.141 would be shown as 3.14.
"""
function mermaid(T::WFST{N,IL,OL,W}; rankdir::String="LR",
                                     keepone ::Bool=true,
                                     keepsame::Bool=true,
                                     digits::Int=2) where {N,IL,OL,W}
    I = T.starts
    F = T.finals
    Q = T.states
    s = "graph $rankdir\n"

    # dealing with arcs leaving src
    for (src, arcs) ∈ T.states
        # handle starts and finals
        bestart = isstart(T, src)
        befinal = isfinal(T, src)
        label   = one(String)
        ID = string(src)
        s *= " style $(ID) "
        if bestart && befinal
            wi = round(I[src]; digits)
            wf = round(F[src]; digits)
            vi = ifelse(isone(wi) && !keepone, "" , "$wi")
            vf = ifelse(isone(wf) && !keepone, "" , "$wf")
            label *= "$ID/$vi,$vf"
            s *= "fill:#F5F5F5,stroke:#FF7F00,stroke-width:4px\n" # 255 165 0, orange1
        elseif bestart
            w = trunc(I[src]; digits)
            v = ifelse(isone(w) && !keepone, "" , "/$w")
            label *= "$ID$v"
            s *= "fill:#F5F5F5,stroke:#8A2BE2,stroke-width:2px\n" # 232 232 232, gray91
        elseif befinal
            w = trunc(F[src]; digits)
            v = ifelse(isone(w) && !keepone, "" , "/$w")
            label *= "$ID$v"
            s *= "fill:#F5F5F5,stroke:#CD3333,stroke-width:2px\n" # 105 105 105, grey41
        else
            label *= "$ID"
            s *= "fill:#F5F5F5,stroke:#000000\n"
        end

        for arc ∈ arcs
            dst = string(next(arc))
            i = ilabel(arc)
            o = olabel(arc)
            w = weight(arc)
            w = trunc(w; digits)
            is = ifelse(isone(i), "ε", string(i))
            os = ifelse(isone(o), "ε", string(o))
            ws = ifelse(isone(w) && !keepone, "" , "$w")
            e = ifelse(is==os && !keepsame, "$is$ws", "$is:$os/$ws")
            s = "$s   $ID([$label]) -- \"$e\" --> $dst([$dst])\n"
        end
    end
    return s * "%%{init:{'flowchart':{'curve':'basis'}}}%%\n"
end
