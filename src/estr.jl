@inline hasblank(e::String) = occursin(" ", e)

"""
    Estr(e::String)

Through the instances of Estr, three kinds of symbolic computation, `+`, `*` and `inv` are provided.
`Estr` is intended for checking the weight expression when appling WFST operations e.g. `compose`, `determine` 
`rmeps` etc with no need of numerical calculations.

# Example
```julia
begin
    @estr a b c d e f g
    T = WFST{Int, String, String, Estr}();
    addstart(T, 0)
    addeps(T, 0,1, b)
    addeps(T, 1,2, c)
    addarc(T, 2,3, "X","X", d)
    addarc(T, 2,4, "Y","Y", e)
    addarc(T, 1,2, "Z","Z", g)
    addeps(T, 2,5, f)
    addfinal(T, 3)
    addfinal(T, 4)
    addfinal(T, 5)
    R = rmeps(T)
    D = connect(R)
    open("./fst.dot", "w") do io
        write(io, wfsts2dot([T,R,D],["T", "R", "D"], keepsame=false));
    end
end
```
"""
struct Estr
    e :: String
    function Estr(e::String)
        isone(e) && error("input is empty :/")
        hasblank(e) && error("input shall not contain any blank :/")
        new(e)
    end
end


function Base.:+(x::Estr, y::Estr)
    iszero(x) && return y
    iszero(y) && return x
    return Estr("($(x.e)+$(y.e))")
end


function Base.:*(x::Estr, y::Estr)
    isone(x) && return y
    isone(y) && return x
    return Estr("$(x.e)*$(y.e)")
end


function Base.show(io::IO, ::MIME"text/plain", x::Estr)
    show(io, x.e)
end


Base.zero(::Type{Estr}) = Estr("0")
Base.one(::Type{Estr})  = Estr("1")
Base.zero(x::Estr) = Estr("0")
Base.one(x::Estr)  = Estr("1")

Base.:(==)(x::Estr, y::Estr)   = isequal(x.e, y.e)
Base.isequal(x::Estr, y::Estr) = isequal(x.e, y.e)

Base.iszero(x::Estr) = Estr("0")==x
Base.isone(x::Estr)  = Estr("1")==x


function Base.inv(x::Estr)
    s = x.e
    if length(s)==1
        return Estr("$(s)⁻¹")
    end
    if s[1] == '('
        return Estr("$(s)⁻¹")
    end
    return Estr("inv($s)")
end


Base.trunc(x::Estr;digits) = x.e
Base.string(x::Estr) = x.e


"""
    @estr x1 x2 x3 ...

generate Estr instances x1 x2 x3 ...

# Example
    julia> @estr x y z;
    julia> x + y * z
    "(x+y*z)"
"""
macro estr(xs...)
    ys = Expr(:block)
    for x ∈ xs
        e = Estr(string(x))
        push!(ys.args, Expr(:(=), x, e))
    end
    esc(ys);
end

