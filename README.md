# Dimers

[![Build Status](https://travis-ci.org/sswatson/Dimers.jl.svg?branch=master)](https://travis-ci.org/sswatson/Dimers.jl)

`Dimers` is a package for simulating the
[dimer model](http://en.wikipedia.org/wiki/Domino_tiling) on a 2D
rectangular grid. It also provides support for loop erased random walks and
Wilson's algorithm on an arbitrary graph.

```julia
showgraphics(draw_graph(dimer_sample(20)))
```

We can also compute the height function associated with the dimer sample:

```julia
dimer_height(dimer_sample(20))
```

`Wilson` takes a graph as its first argument and an array of `true`/`false`
values specifying the roots. 

```julia
showgraphics(draw_graph(Wilson(G,[[true];[false for i=1:length(G.vertices)-1]])))
```

`LERW(G,v0,roots)` samples a loop-erased random walk on the graph `G`
starting from the vertex whose index in `G.vertices` is `v0` stopped upon
hitting one of the vertices `v` for which `roots[v]` is `true`. 

```julia
n = 100
G = Graphs.adjlist((Int64,Int64),is_directed=false)

for i=1:n
    for j=1:n
        Graphs.add_vertex!(G,(i,j))
    end
end

roots = Bool[v[1] == 1 || v[1] == n || v[2] == 1 || v[2] == n  for v in
G.vertices];

v0 = find(x->x==(div(n,2),div(n,2)),G.vertices)[1]

lerw = LERW(grid_graph(n),v0,roots)
for i=1:length(lerw)-1
    add_edge!(G,lerw[i],lerw[i+1])
end

showgraphics(draw_graph(G))
```
