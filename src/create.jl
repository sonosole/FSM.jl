"""
    linearwfsa(text::String; dlm::String=" ", wtype::DataType=Float64)
return a linear **WFSA** given by the labels in `text`

# Example
```julia
julia> linearwfsa("hello world")
WFST{Float64}
 starts:
  0/1.0
 finals:
  2/1.0
 transitions:
  0→1 hello:hello/1.0
  1→2 world:world/1.0
```
"""
function linearwfsa(text::String; dlm::String=" ", wtype::DataType=Float64)
    A = WFST{Int, String, String, wtype}()
    t = string.(split(text, dlm))
    for (i, x) ∈ enumerate(t)
        addarc(A, i-1, i, x, x)
    end
    addstart(A, 0)
    addfinal(A, length(t))
    return A
end


"""
    txt2wfsa(txtfile::String; wtype::DataType=Float64)
return a linear **WFSA** described in each line of `txtfile`

# Example
suppose we have a txt file with context:\n
    ni hao
    hi there
then we have the example 
```julia
julia> txt2wfsa("./dict.txt")
WFST{Float64}
 starts:
  0/1.0
 finals:
  4/1.0
  2/1.0
 transitions:
  0→1 ni:ni/1.0
  0→3 hi:hi/1.0
  3→4 there:there/1.0
  1→2 hao:hao/1.0
```
"""
function txt2wfsa(txtfile::String; wtype::DataType=Float64)
    T = WFST{Int, String, String, wtype}()
    b = 0
    for line ∈ eachline(txtfile)
        words = string.(split(line, " "))
        count = length(words)

        addarc(T, 0, b+1, first(words), first(words))
        for i = 2 : count
            j = i + b
            addarc(T, j-1, j, words[i], words[i])
        end
        b += count
        addfinal(T, b)
    end
    addstart(T, 0)
    return T
end


"""
    dictwfst(text; wtype::DataType=Float64)
create a WFST for a dict of one word, where `wtype` is the data type of WFST's weight. 

# Example
```julia
julia> dictwfst("it ih t")
WFST{Float64}
 starts:
  0/1.0
 finals:
  2/1.0
 transitions:
  0→1 ih:it/1.0
  1→2 t:ϵ/1.0
```
"""
function dictwfst(text::String; wtype::DataType=Float64)
    wordspells = string.(split(text, " "))
    n = length(wordspells)
    word = wordspells[1]
    spells = wordspells[2:n]
    nspells = n - 1

    T = WFST{Int, String, String, wtype}()
    ε = ""
    
    addarc(T, 0, 1, first(spells), word)
    for i = 2:nspells
        addarc(T, i-1, i, spells[i], ε)
    end

    addstart(T, 0)
    addfinal(T, nspells)
    return T
end


"""
    dictswfst(txtfile; wtype::DataType=Float64)
Create a dict WFST from `txtfile` recording mutiple words and corresponding pronunciations. 

# Example
suppose we have a dict `txtfile` with:\n
    it ih t
    stop s t aa p
then we have the example code snippet 

```julia
julia> dictswfst("./dict.txt")
WFST{Float64}
 starts:
  0/1.0
 finals:
  6/1.0
  2/1.0
 transitions:
  0→1 ih:it/1.0
  0→3 s:stop/1.0
  4→5 aa:ϵ/1.0
  5→6 p:ϵ/1.0
  3→4 t:ϵ/1.0
  1→2 t:ϵ/1.0
```
"""
function dictswfst(txtfile::String; wtype::DataType=Float64)
    T = WFST{Int, String, String, wtype}()
    ε = ""
    b = 0
    for line ∈ eachline(txtfile)
        wordspells = string.(split(line, " "))
        n = length(wordspells)
        word = wordspells[1]
        spells = wordspells[2:n]
        nspells = n - 1

        addarc(T, 0, b+1, first(spells), word)
        for i = 2:nspells
            j = i + b
            addarc(T, j-1, j, spells[i], ε)
        end

        b += nspells
        addfinal(T, b)
    end
    addstart(T, 0)
    return T
end

