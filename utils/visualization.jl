function plotSolution(graph, x, enum_cycles, cycles_vertic,membership_node, membership_edge)
    n_node = nv(graph)
    n_edge = ne(graph)

    real_cycles = x.==1
    
    nodelabel = 1:n_node
    edgelabel = 1:n_edge
    
    # permet de définir une taille pour les noeuds
    nodesize = [Graphs.outdegree(kep_graph, v) for v in Graphs.vertices(kep_graph)] # nodesize=nodesize


    nodecolor = [colorant"lightgrey", colorant"orange"]
    nodefillc = nodecolor[membership_node]

    edgestrokec = nodecolor[membership_edge]


    return gplot(kep_graph,
        nodelabel=nodelabel,
        nodelabeldist=1.5,
        nodelabelangleoffset=π/4,
        nodefillc=nodefillc,
        edgestrokec=edgestrokec,
        layout=circular_layout)
end


;