module Dimers

export dimers, Wilson, LERW, GridGraph, rotate, flatten, midpoint, dimer_sample, draw_graph, dimer_height

import Graphs, Graphics2D

import Base.getindex

getindex{V}(S::Set{V},k::Integer) = collect(keys(S.dict))[k]

function LERW{V,E}(Γ::Graphs.AbstractGraph{V,E},startingvertex::Int64,roots::Array{Bool,1})
    X = [startingvertex]
    maxiter = 10^7
    cntr = 1
    while ~(roots[X[end]])
        neighbors = Γ.adjlist[X[end]]
        push!(X,Graphs.vertex_index(neighbors[rand(1:length(neighbors))],Γ))
        i = 1
        while X[i] != X[end]
            i += 1
        end
        X = X[1:i]
        cntr += 1
        if cntr >= maxiter
            error("Maximum iterations hit in LERW")
        end
    end
    vert = Graphs.vertices(Γ)
    return vert[X]
end

function Wilson{V,E}(Γ::Graphs.AbstractGraph{V,E},roots::Array{Bool,1})
    if length(Graphs.connected_components(Γ)) > 1
        error("Graph not connected")
    end
    maxiter = length(Graphs.vertices(Γ))
    cntr = 1
    UST = Graphs.adjlist(typeof(Graphs.vertices(Γ)[1]))
    for v in Graphs.vertices(Γ)
        Graphs.add_vertex!(UST,v)
    end
    discovered = roots
    while ~all(discovered)
        i = 1
        while discovered[i]
            i += 1
        end
        lerw = LERW(Γ,i,discovered)
        for k=1:length(lerw)-1
            Graphs.add_edge!(UST,lerw[k+1],lerw[k])
            discovered[Graphs.vertex_index(lerw[k],Γ)] = true
        end
        cntr += 1
        if cntr >= maxiter
            error("Something's gone wrong with Wilson's algorithm")
            break
        end
    end
    return UST
end

rotate(segment::((Int64,Int64),(Int64,Int64))) = 
(div(segment[1][1] + segment[2][1] - segment[1][2] + segment[2][2],2), 
div(segment[1][1] - segment[2][1] + segment[1][2] + segment[2][2],2)),
(div(segment[1][1] + segment[2][1] + segment[1][2] - segment[2][2],2), 
div(-segment[1][1] + segment[2][1] + segment[1][2] + segment[2][2],2))

flatten(e::((Int64,Int64),(Int64,Int64))) = vcat(map(x->vcat(x...),vcat(e...))...)

midpoint(point1::(Int64,Int64),point2::(Int64,Int64)) = (div(point1[1] + point2[1],2),div(point1[2] + point2[2],2))

function GridGraph(n::Int64)

    Γ = Graphs.adjlist((Int64,Int64),is_directed=false)

    for i=1:n
        for j=1:n
            Graphs.add_vertex!(Γ,(i,j))
        end
    end

    gridedges = ((Int64,Int64),(Int64,Int64))[]

    for i=1:n
        for j=1:n
            if i < n
                Graphs.add_edge!(Γ,(i,j),(i+1,j))
            end
            if j < n
                Graphs.add_edge!(Γ,(i,j),(i,j+1))
            end
        end
    end
    
    return Γ 
end

function dimer_sample(m::Int64,n::Int64)
    
    function add_edge_and_continue(vertex::(Int64,Int64),prev_vertex::(Int64,Int64))
        Graphs.add_edge!(dualtree_ordered,prev_vertex,vertex)
        if (vertex[1]+2,vertex[2]) != prev_vertex && (vertex[1]+2,vertex[2]) in Graphs.out_neighbors(vertex,dualtree)
            add_edge_and_continue((vertex[1]+2,vertex[2]),vertex)
        end
        if (vertex[1]-2,vertex[2]) != prev_vertex && (vertex[1]-2,vertex[2]) in Graphs.out_neighbors(vertex,dualtree)
            add_edge_and_continue((vertex[1]-2,vertex[2]),vertex)
        end
        if (vertex[1],vertex[2]+2) != prev_vertex && (vertex[1],vertex[2]+2) in Graphs.out_neighbors(vertex,dualtree)
            add_edge_and_continue((vertex[1],vertex[2]+2),vertex)
        end
        if (vertex[1],vertex[2]-2) != prev_vertex && (vertex[1],vertex[2]-2) in Graphs.out_neighbors(vertex,dualtree)
            add_edge_and_continue((vertex[1],vertex[2]-2),vertex)
        end
    end

    Γ = Graphs.adjlist((Int64,Int64),is_directed=false)

    for i=1:m+1
        for j=1:n+1
            if (i,j) != (m+1,n+1)
                Graphs.add_vertex!(Γ,(2i-1,2j-1))
            end
        end
    end

    primaledges = ((Int64,Int64),(Int64,Int64))[]

    for i=1:m
        for j=1:n
            push!(primaledges,((2i-1,2j-1),(2i-1+2,2j-1)))
            push!(primaledges,((2i-1,2j-1),(2i-1,2j-1+2)))
        end
    end

    for e in primaledges
        Graphs.add_edge!(Γ,e...) 
    end

    dualtree = Graphs.adjlist((Int64,Int64),is_directed=false)

    for i=0:m
        for j=0:n
            if (i,j) != (0,0)
                Graphs.add_vertex!(dualtree,(2i,2j))
            end
        end
    end

    roots = [false for i=1:length(Graphs.vertices(Γ))]
    dualroots = [false for i=1:length(Graphs.vertices(dualtree))]

    for i=1:m
        roots[Graphs.vertex_index((2i-1,2n+1),Γ)] = true
        dualroots[Graphs.vertex_index((2i,0),dualtree)] = true
    end

    for j=1:n
        roots[Graphs.vertex_index((2m+1,2j-1),Γ)] = true
        dualroots[Graphs.vertex_index((0,2j),dualtree)] = true
    end

    UST = Wilson(Γ,roots)

    for e in primaledges
        if ~(e[2] in Graphs.out_neighbors(e[1],UST) || e[1] in Graphs.out_neighbors(e[2],UST))
            newedge = rotate(e)
            Graphs.add_edge!(dualtree,newedge...)
        end
    end

    dualtree_ordered = Graphs.inclist((Int64,Int64),is_directed=true)

    for v in Graphs.vertices(dualtree)
        Graphs.add_vertex!(dualtree_ordered,v) 
    end

    for k=2:2:2n
        if (2,k) in Graphs.out_neighbors((0,k),dualtree) || (0,k) in Graphs.out_neighbors((2,k),dualtree)
            add_edge_and_continue((2,k),(0,k))
        end
        if (k,0) in Graphs.out_neighbors((k,2),dualtree) || (k,2) in Graphs.out_neighbors((k,0),dualtree)
            add_edge_and_continue((k,2),(k,0))
        end    
    end

    dimergraph = Graphs.adjlist(typeof(Graphs.vertices(Γ)[1]))
    dimer_vertices = (Int64,Int64)[]
    dimer_edges = ((Int64,Int64),(Int64,Int64))[]

    for i=1:2m
        for j=1:2n
            push!(dimer_vertices,(i,j))
            if i < 2m
                push!(dimer_edges,((i,j),(i+1,j)))
            end
            if j < 2n
                push!(dimer_edges,((i,j),(i,j+1)))
            end
            Graphs.add_vertex!(dimergraph,(i,j))
        end
    end

    for v in Graphs.vertices(UST)
        for w in Graphs.out_neighbors(v,UST)
            Graphs.add_edge!(dimergraph,midpoint(v,w),w)
        end
    end

    for v in Graphs.vertices(dualtree_ordered)
        for w in Graphs.out_neighbors(v,dualtree_ordered)
            Graphs.add_edge!(dimergraph,midpoint(v,w),w)
        end
    end

    return dimergraph

end

dimer_sample(n::Integer) = dimer_sample(n,n)

function draw_graph{V,E}(Γ::Graphs.AbstractGraph{V,E})
    
    all_points = Graphics2D.GraphicElement[]
    all_edges = Graphics2D.GraphicElement[]
    
    for v in Graphs.vertices(Γ)
        push!(all_points,Graphics2D.Point(v...))
        for w in Graphs.out_neighbors(v,Γ)
            push!(all_edges,Graphics2D.Line([v[1] v[2]; w[1] w[2]];rs=1.0))
        end
    end
    
    return [all_points; all_edges]
end

function dimer_height{V,E}(dimer_graph::Graphs.AbstractGraph{V,E})
    
    all_edges = ((Int64,Int64),(Int64,Int64))[]
    
    m = div(maximum(map(x->x[1],Graphs.vertices(dimer_graph))),2)
    n = div(maximum(map(x->x[2],Graphs.vertices(dimer_graph))),2)

    for i=1:2m
        for j=1:2n
            if i < 2m
                push!(all_edges,((i,j),(i+1,j)))
            end
            if j < 2n
                push!(all_edges,((i,j),(i,j+1)))
            end
        end
    end

    down_edges = zeros(Bool,2m+1,2n+1)

    for i=1:2m-1
        for j=1:2n
            if (i,j) in Graphs.out_neighbors((i+1,j),dimer_graph) || (i+1,j) in Graphs.out_neighbors((i,j),dimer_graph)
                down_edges[i,j] = true
            end        
        end
    end

    h = zeros(Int64,2n+1,2n+1)

    for j=2:2:2n+1
        h[1,j] = 1
    end

    for i=2:2n+1
        h[i,1] = isodd(i) ? 0 : 1
        for j=2:2n+1
            if down_edges[i-1,j-1]
                h[i,j] = h[i-1,j] + (-1)^(i+j)
            else
                h[i,j] = h[i,j-1] + (-1)^(i+j+1)
            end
        end
    end
    
    return h
end

end # module
