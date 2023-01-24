
"""
addSpecificEdges

Add specific edges

# Parameters
* `kep_graph::MetaDiGraph` : graph describing the pairs and compatibilities
* `list_edges` : edges to add
* `fail` : failure rate to affect
* `weight` : weight to affect
"""
function addSpecificEdges(kep_graph, list_edges, fail, weight)
    n_row = size(list_edges)[1]
    print(n_row)
    for i in 1:1:n_row
        add_edge!(kep_graph, list_edges[i, 1], list_edges[i, 2], :weight, weight)
        set_prop!(kep_graph, Edge(list_edges[i, 1], list_edges[i, 2]), :failure, fail)
    end
    return kep_graph
end;

"""

"""
function addRandomEdge(kep_graph, fail, weight)
    not_stop = true
    nb_nodes = ne(kep_graph)
    while(not_stop)
        s = trunc(Int, rand(1:nb_nodes))
        d = trunc(Int, rand(1:nb_nodes))
        cond = ((s != d) && (~has_edge(kep_graph, s, d)))
        if cond
            add_edge!(kep_graph, s, d, :weight, weight)
            set_prop!(kep_graph, Edge(s, d), :failure, fail)
            not_stop = false
        end
    end
end


"""
    addRandomEdges

Add random edges

# Parameters
* `kep_graph::MetaDiGraph` : graph describing the pairs and compatibilities
* `ratio` : ratio 
* `fail` : failure rate to affect
* `weight` : weight to affect
"""
function addRandomEdges(kep_graph, ratio, fail, weight)
    nb_nodes = nv(kep_graph)
    nb_edges = ne(kep_graph)
    nb_it = floor((((nb_nodes+1) * nb_nodes)/2) * ratio)
    for i in 1:1:nb_it
        addRandomEdge(kep_graph, fail, weight)
    end
    
    return kep_graph
end
; 

function addRandomCycles(kep_graph, nb_cycles, fail, weight, size_cluster=3)
    curr_nb_cycles = length(simplecycles_limited_length(kep_graph, size_cluster, 10^6))
    #print(curr_nb_cycles,  nb_cycles)
    while length(simplecycles_limited_length(kep_graph, size_cluster, 10^6)) < curr_nb_cycles + nb_cycles
        #println(length(simplecycles_limited_length(kep_graph, size_cluster, 10^6)))
        addRandomEdge(kep_graph, fail, weight)
    end
    return kep_graph
end
;


function removeUselessNodes(kep_graph, size_cluster=3)
    nb_nodes = nv(kep_graph)
    matrix_nodes_temp = Int64[j for i in 1:nb_nodes, j in 1:nb_nodes]
    matrix_nodes = Int64[j for i in 1:nb_nodes, j in 1:nb_nodes]
    for i in 1:1:nb_nodes
        matrix_nodes_temp[i, :] = dijkstra_shortest_paths(kep_graph, i).dists
    end

    for i in 1:1:nb_nodes
        for j in 1:1:nb_nodes
            if matrix_nodes_temp[i, j]>size_cluster || matrix_nodes_temp[j, i]>size_cluster || matrix_nodes_temp[i, j] + matrix_nodes_temp[j, i]>size_cluster
                matrix_nodes[i, j] = 0
            else
                matrix_nodes[i, j] = matrix_nodes_temp[i, j] + matrix_nodes_temp[i, j] 
            end
            
        end
    end
    matrix_nodes = sum(matrix_nodes, dims=2)
    return matrix_nodes
end