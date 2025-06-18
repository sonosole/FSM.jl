# module FSM

include("type-arc.jl")
export Arc
export swapio!, iseps, isieps, isoeps, ismatched
export ilabel, olabel, weight, next

include("type-fst.jl")
export WFST
export isstart, isfinal
export haseps, hasieps, hasoeps, hasioeps
export iswfsa
export starts, finals
export ilabelsof, olabelsof, iolabelsof
export nstates, nstarts, nfinals, narcs, neps, nieps, noeps, nioeps

include("op-basic.jl")
export addstart, addfinal, addeps, addarc, addpath
export rmstart, rmfinal, rmstate

include("order.jl")
export keyorder, reorder


###########################
# normalization operations
###########################
include("norm-start-final.jl")
include("norm-push-weight.jl")
include("norm-sort-arcs.jl")
include("norm-sort-topo.jl")
export onestart, onefinal
export potential, pushweight, reweight, push
export arcsort!
export indegree, toposort, levelsort

###########################
# rational operations
###########################
include("rational-concat.jl")
include("rational-star-plus.jl")
include("rational-union.jl")
export star, star!, plus, plus!

###########################
# unary operations
###########################
include("unary-invert.jl")
include("unary-project.jl")
include("unary-reverse.jl")
export swapio, swapio!, invert, invert!
export iiproj, iiproj!, ooproj, ooproj!

###########################
# binary operations
###########################
include("binary-compose.jl")
include("binary-intersect.jl")
export compose, _compose

###########################
# search operations
###########################
include("shortest-dist-eps.jl")
include("shortest-dist-symbol.jl")
include("shortest-dist.jl")
include("shortest-path.jl")
export dists, dist, shortestdist
export epsdists, ϵdists
export symdists

###########################
# optimization operations
###########################
include("optimal-connect.jl")
include("optimal-detwfsa.jl")
include("optimal-detwfst.jl")
include("optimal-minimize.jl")
include("optimal-rmepsilon.jl")
export connect
export adet, detwfsa, isdet, issequential, isfunctional
export tdet, detwfst, det
export rdrdfsa, rdrd
export epsqv, ϵindegree, ϵtoposort, until_not_eps_states
export rmeps, rmϵ, rm_acyclic_eps

include("op-combination.jl")
export ideti

###########################
# pretty showing operations
###########################
include("io-mermaid.jl")
include("io-graphviz.jl")
export mermaid
export dot, dots, savedot, savedots, savewfst, leveldot, saveleveldot, linedot

# end # module FSM
