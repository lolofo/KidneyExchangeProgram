function getCluster(kep_graph, sol)
    nodelabel = [v for v in Graphs.vertices(kep_graph)]
    stacklabel = []
    cluster = Dict()
    k = 1
    for i in nodelabel
        if i ∉ stacklabel
            cluster[k] = [i]
            for j in nodelabel[(i+1):length(nodelabel)]
                if sol[i, j] == 1 
                    append!(stacklabel, j)
                    append!(cluster[k], j)
                end
            end
        end
        k += 1
    end
    return cluster
end;


function getClusterUsefull(cluster)
    clusterUsefull = Dict()
    for (key, value) in cluster
        if length(value)>1
            clusterUsefull[key] = value
        end
    end
    return clusterUsefull
end;
