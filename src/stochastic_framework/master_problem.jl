
using GLPK

""" 
    masterProblem

# Parameters
* `C` : index of the cycles
* `vertic_cycles` : a dictionnary, at the key i of this list, give the a list of the index of the cycles
                    which involve the node i
* `U` : the utilities of each cycle in the graph.
"""
function masterProblem(C, vertic_cycles ,U)

    model = Model(GLPK.Optimizer)

    @variable(model, x[i in C], Bin)

    @objective(model, Min, sum(- x[c]U[c] for c in C))

    for (v, C_v) in vertic_cycles
        @constraint(model, sum(x[c] for c in C_v)<=1)
    end
    
    return(Dict("model" => model))

end
;

"""
    addThetaVar

# Parameters
* `model` : a models which comes from the master problem type
* `K` : the number of scenarios
"""
function addThetaVar(model, K)

    @variable(model, theta[i in 1:1:K]>=0)

    return(Dict("model" => model))
end
;

"""
    updateMasterProblem

Update the master problem to have the classical recourse modelisation

# Parameters
* `model` : a master problem type
* `duals` : the duals solutions associated with each scenarios (of shape |C| x S)
* `T_array` : the different scenarios (of shape |C| x S)
* `C` : our index cycles
* `U` : the utility of each cycle in our cycles (of shape |C|)
"""
function updateMasterProblem(model, duals, T_array, C, U)

    S = size[T_array][2] # get the number of scenarios

    @objective(model, Min, sum(- x[c]U[c] for c in C) + (1/S) * sum(theta[s] for s in 1:1:S))

    for s in 1:1:S
        @constraint(model, theta[s] >= sum(T_array[c, s] * x[c] * duals[c, s] for c in C))
    end

    return(Dict("model" => model))
end
;



"""
    updateCVARMasterProblem

Update the master problem to have an risk adverse modelisation based on the CVAR criterion

# Parameters
* `model` : a master problem type
* `alpha` : the hazard rate for the CVaR modelisation
* `duals` : the duals solutions associated with each scenarios (of shape |C| x S)
* `T_array` : the different scenarios (of shape |C| x S)
* `C` : our index cycles
* `U` : the utility of each cycle in our cycles (of shape |C|)
"""
function updateCVARMasterProblem(model, alpha, duals, T_array, C, U)
    # TODO
    return Nothing
end
;


################################################################################################
################################################################################################
################################################################################################



"""
In this part we propose a new insiste on the kindey exchange programe based on cluster
"""

""" 
    masterClusterProblem

This function will initiate the master problem.
The master problem try to find the best clusters for our problem.

# Parameters
* `kep_graph` : the kep graph
* `ClusterSize` (int) : the maximal size of the clusters

# Return 
* this function will return the master problem coded with its constraints but with no objective
"""
function masterClusterProblem(kep_graph, ClusterSize)

    V = vertices(kep_graph) # the vertices of our graph

    
    model = Model(GLPK.Optimizer)
    @variable(model, x[i = V, j = V], Bin)
    
    # lets define the clustering constraints
    for i in V
        for j in V
            # réflexivity
            @constraint(model, x[i, j] == x[j, i])
            for k in V
                # transitivity
                @constraint(model, x[i, j]+x[j, k]-1 <= x[i, k])
            end
        end
        # size of the clusters
        @constraint(model, sum(x[i, j] for j in V) <= ClusterSize)
    end
    return(Dict("model" => model))
end
;