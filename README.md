# Dimers

`Dimers` is a package for simulating the
[dimer model](http://en.wikipedia.org/wiki/Domino_tiling) on a 2D
rectangular grid, using
[an algorithm of Kenyon, Propp, and Wilson](http://arxiv.org/pdf/math/9903025v2.pdf). `Dimers`
also provides support for loop erased random walks and Wilson's algorithm
on an arbitrary graph.

```julia
showgraphics(draw_graph(dimer_sample(20)))
```

![Dimer sample](https://github.com/sswatson/Dimers.jl/blob/master/images/dimersample.png)

We can also compute the height function associated with the dimer sample:

```julia
dimer_height(dimer_sample(20))
```

```
11x11 Array{Int64,2}:
  0   1   0   1   0  1   0   1   0   1   0
 -1  -2  -1  -2  -1  2  -1  -2  -1  -2  -1
  0  -3  -4  -3   0  1   0   1   0  -3   0
 -1  -2  -1  -2  -1  2  -1  -2  -1  -2  -1
  0  -3   0   1   0  1   0   1   0   1   0
 -1  -2  -1   2  -1  2  -1  -2  -1   2  -1
  0   1   0   1   0  1   0   1   0   1   0
 -1  -2  -1   2  -1  2   3   2  -1  -2  -1
  0   1   0   1   0  1   0   1   0   1   0
 -1   2   3   2   3  2   3   2   3   2  -1
  0   1   0   1   0  1   0   1   0   1   0
```

`Wilson` takes a graph as its first argument and an array of `true`/`false`
values specifying the roots. 

```julia
showgraphics(draw_graph(Wilson(G,[[true];[false for i=1:length(G.vertices)-1]])))
```

![UST sample](https://github.com/sswatson/Dimers.jl/blob/master/images/USTsample.png)

`LERW(G,v0,roots)` samples a loop-erased random walk on the graph `G`
starting from the vertex whose index in `G.vertices` is `v0` stopped upon
hitting one of the vertices `v` for which `roots[v]` is `true`. 

```julia
import Graphs
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

![Loop-erased random walk sample](https://github.com/sswatson/Dimers.jl/blob/master/images/lerwsample.png)

[![Build Status](https://travis-ci.org/sswatson/Dimers.jl.svg?branch=master)](https://travis-ci.org/sswatson/Dimers.jl)
