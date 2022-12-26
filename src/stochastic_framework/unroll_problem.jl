using Random, MetaGraphs, SimpleWeightedGraphs, Graphs, JuMP, DelimitedFiles, Distributions, Plots, GraphPlot
""" 
    unRollPRoblem

This function has the objective to solve the problem of the kidney exchange program buy unrolling all the scenarios.

# Parameters
* `kep_graph` : the kep graph
* `ClusterSize` (int) : the maximal size of the clusters
* `C` : the index of the cycles
* `cycles` : the array of the cycles of length <= k
* `U` : the vector of the utility of the cycles
* `ksi_s` : a 3D tensor which gives us the scenario, if we want ksi_{i,j}^k the value of ksi_{i,j} for the scenario k, is
ksi_s[i,j,k]

# Return
The unroll model, ready to be optimize.

"""
function unrollClusterProblem(kep_graph, ClusterSize, C, cycles, U, ksi_s, vertic_cycles)

    V = vertices(kep_graph) # the vertices of our graph
    K = 1:1:size(ksi_s)[3] # the different scenarios
    
    model = Model(GLPK.Optimizer)

    # the cluster variables
    @variable(model, x[i = V, j = V], Bin)
    @variable(model, y[c = C, k = K], Bin)
    @objective(model, Max, (1/length(K)) * sum(sum(y[c,k]*U[c] for c in C) for k in K))
    
    # lets define the clustering constraints
    for i in V
        for j in V
            # reflexivity
            @constraint(model, x[i, j] == x[j, i])
            for k in V
                # transitivity
                @constraint(model, x[i, j]+x[j, k]-1 <= x[i, k])
            end
        end
        # size of the clusters
        @constraint(model, sum(x[i, j] for j in V) <= ClusterSize)
    end

    # unroll all the constraints for all the scenarios
    for scenar in K
        for c in C
            current_cycle = cycles[c]
            for k in 1:1:length(current_cycle)
                i = current_cycle[k]
                j = Nothing # for initialisation purposes
                if k==length(current_cycle)
                    # if we are at the end the following node is the first one
                    j = current_cycle[1] 
                else
                    j = current_cycle[k+1] # following node in the cycle
                end
                # we can choose a cycle iif the edge exists and the cross test is successful
                @constraint(model, y[c, scenar] <= x[i, j]*ksi_s[i, j, scenar])
            end 
        end
        
        for (v, C_v) in vertic_cycles
            # each node must be at most in one cycle
            @constraint(model, sum(y[c, scenar] for c in C_v)<=1)
        end
    end
    return(Dict("model" => model))
end
;