
"""
getCluster

This function allow us to exctract clusters, cluster with one element are considered
# Parameters

*`graph` : the kep graph
*`sol` : The solution of the master problem
"""

function getCluster(graph, sol)
    nodelabel = [v for v in Graphs.vertices(graph)]
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

"""
getClusterUsefull

This function allow us to exctract clusters, cluster with less than 2 elements aren't considered
# Parameters

*`cluster` : solution given by getCluster(graph, sol))
"""
function getClusterUsefull(cluster)
    clusterUsefull = Dict()
    for (key, value) in cluster
        if length(value)>1
            clusterUsefull[key] = value
        end
    end
    return clusterUsefull
end;
