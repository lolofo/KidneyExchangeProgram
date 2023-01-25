
using GLPK;
using HiGHS;


### the clustering method ###
#############################

""" 
    masterClusterProblem

This function will initiate the master problem.
The master problem here will try to solve the mean value problem.
Thanks to this function we will have access to the x_{EV}, convenient to calculate the VSS

# Parameters
* `kep_graph` : the kep graph
* `ClusterSize` (int) : the maximal size of the clusters
* `C` : the index of the cycles
* `cycles` : the array of the cycles of length <= k
* `U` : the vector of the utility of the cycles

# Return 
* this function will return the master problem coded with its constraints
"""
function masterClusterProblem(kep_graph, ClusterSize, C, cycles, U, vertic_cycles)

    V = vertices(kep_graph) # the vertices of our graph
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # the cluster variables
    @variable(model, x[i = V, j = V], Bin)

    # the mean value problem variable.
    @variable(model, 1 >= z[c in C] >=0)

    @objective(model, Max, sum(z[c]U[c] for c in C))
    
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

    # let's define the constraint for the mean value problem
    for c in C
        current_cycle = cycles[c]
        for k in 1:1:length(current_cycle)
            i = current_cycle[k]
            j = Nothing
            if k==length(current_cycle)
                j = current_cycle[1] 
            else
                j = current_cycle[k+1]
            end
            @constraint(model, z[c] <= x[i, j]*(1 - get_prop(kep_graph, Edge((i,j)), :failure)))
        end 
    end

    for (v, C_v) in vertic_cycles
        # each node must be at most in one cycle
        @constraint(model, sum(z[c] for c in C_v)<=1)
    end

    return(Dict("model" => model))
end
;



### the updates of the master problem for the ###
#################################################

"""
    addThetaCluster

This function is part of the l-shape method. It will modify the objective of the cluster problem.

# Parameters
*`model` : the master problem
* `nb_scenar` : the number of scenarios
"""
function addThetaCluster(model, nb_scenar, C, U)

    V = 1:1:size(model[:x])[1]
    # add the new var
    @variable(model, theta[k in 1:1:nb_scenar])

    # update the objective
    # sum(model[:z][c]U[c] for c in C) + 
    @objective(model, Max, (1/nb_scenar)*sum(model[:theta][k] for k in 1:1:nb_scenar))
end;


"""
    updateCluster

This function will add the L-type constraints, to the master masterProblem

# Parameters
* `kep_graph` : the graph of the kidney exchange problem
* `model` : the master problem
* `C`: the list of the cycle index
* `ksi_k` : the scenario, it's an array of shape (|V|, |V|)
* dual : a dictionnary containing the dual solution corresponding to the scenario ksi_k
* `k` : the index of the scenario
"""
function updateCluster(kep_graph, model, C, ksi_k, dual, k)
    V = nv(kep_graph)

    lambda = dual["dual_lambda"]
    mu = dual["dual_mu"]
    delta = dual["dual_delta"]

    @constraint(model, model[:theta][k] <= sum(sum(sum(lambda[i,j,c]*model[:x][i,j]*ksi_k[i,j] for i in 1:1:V) for j in 1:1:V) for c in C) + sum(mu[i] for i in 1:1:V) + sum(delta[c] for c in C))
end
;


# risk averse optimization
##########################

"""
    addCVaRVariables

This function will be usefull for the 

# Parameters
*`model` : the master problem
* `nb_scenar` : the number of scenarios
* `alpha` : float for the risk level
"""
function addCVaRVariables(model, nb_scenar, alpha)

    V = 1:1:size(model[:x])[1]
    # add the new var
    @variable(model, t)
    @variable(model, theta[k in 1:1:nb_scenar])
    @variable(model, pi_var[k in 1:1:nb_scenar]>=0)

    for k in 1:1:nb_scenar
        @constraint(model, model[:pi_var][k] >= -model[:theta][k] - model[:t])
    end

    @objective(model, Min, model[:t] + 1/(1-alpha) * (1/nb_scenar)*sum(model[:pi_var][k] for k in 1:1:nb_scenar))
end
;