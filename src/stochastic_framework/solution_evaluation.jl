using HiGHS;
using GLPK;

include(join(["recourse_problem.jl"], Base.Filesystem.path_separator))
include(join(["..", "..", "utils", "monte_carlo.jl"], Base.Filesystem.path_separator))

"""
    evaluateSolution_ls

This function allow us to evaluate the L-shape solution

# Parameters

*`kep_graph` : the kep graph
*`ξ` : the scenarios for the evaluation
*`x` : l-shaped solution value
*`C` : our index cycles
* `vertic_cycles` : a dictionnary, at the key i of this list, give the a list of the index of the cycles
                    which involve the node i
*`U` : the utility of each cycle in our cycles (of shape |C|)
*`cycles` : the array of the cycles of length <= k
"""
function evaluateSolution_ls(kep_graph, ξ, x, C, vertic_cycles, U, cycles)
    nb_scenar = size(ξ)[3]
    j = 0
    obj = 0
    nb_pers_cluster = 0
    nb_cluster = 0
    nb_pers_transp = 0
    nb_cycle_transp = 0
    nb_pers_by_cycle = length.(cycles[C])
    model = nothing

    for i in 1:1:(nb_scenar)
        if j == 0
            j += 1 
            model = recourseClusterProblem(x, ξ[:, :, i], C, vertic_cycles, U, cycles)
            optimize!(model)
            obj += objective_value(model)
            
            
            res = getClusterUsefull(getCluster(kep_graph, x))
            nb_pers_cluster += sum([length(res[k]) for k in keys(res)])
            nb_cluster += length(res)
            nb_pers_transp += sum(value.(model[:y])[:, 1].*nb_pers_by_cycle)
            nb_cycle_transp += sum(value.(model[:y])[:, 1])
            
        else
            
            modifyRecourseClusterProblem(model, x, C, cycles, ξ[:, :, i])
            optimize!(model)
            obj += objective_value(model)
            
            
            res = getClusterUsefull(getCluster(kep_graph, x))
            nb_pers_cluster += sum([length(res[k]) for k in keys(res)])
            nb_cluster += length(res)
            nb_pers_transp += sum(value.(model[:y])[:, 1].*nb_pers_by_cycle)
            nb_cycle_transp += sum(value.(model[:y])[:, 1])
        end
    end
    return Dict(
        "z_sp" => obj / nb_scenar,
        "nb_pers_cluster" => nb_pers_cluster / nb_scenar,
        "nb_cluster" => nb_cluster / nb_scenar,
        "nb_pers_transp" => nb_pers_transp / nb_scenar,
        "nb_cycle_transp" => nb_cycle_transp / nb_scenar
    )
end
;


"""
    evaluateSolution_ws

This function allow us to evaluate the wait and see problem

# Parameters

*`kep_graph` : the kep graph
*`ξ` : the scenarios for the evaluation
*`C` : our index cycles
* `vertic_cycles` : a dictionnary, at the key i of this list, give the a list of the index of the cycles
                    which involve the node i
*`U` : the utility of each cycle in our cycles (of shape |C|)
*`cycles` : the array of the cycles of length <= k
"""
function evaluateSolution_ws(kep_graph, ξ, C, vertic_cycles, U, cycles, ClusterSize)
    nb_scenar = size(ξ)[3]
    j = 0
    obj = 0
    model = Nothing
    for i in 1:1:(nb_scenar)

        model = clusterProblem_ws(kep_graph, ξ[:, :, i], ClusterSize, C, cycles, U, vertic_cycles)
        optimize!(model)
        obj += objective_value(model)

    end
    return obj / nb_scenar
end
;


"""
    clusterProblem_ws

This function crite the wait and see problem

# Parameters

*`kep_graph` : the kep graph
*`ksi` : one scenario
*`ClusterSize`: cluster size
*`C` : our index cycles
*`cycles` : the array of the cycles of length <= k
*`U` : the utility of each cycle in our cycles (of shape |C|)
* `vertic_cycles` : a dictionnary, at the key i of this list, give the a list of the index of the cycles
                    which involve the node i
"""
function clusterProblem_ws(kep_graph, ksi, ClusterSize, C, cycles, U, vertic_cycles)

    V = vertices(kep_graph) # the vertices of our graph
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # the cluster variables
    @variable(model, x[i = V, j = V], Bin)

    # the mean value problem variable.
    @variable(model, y[c in C] >= 0)

    @objective(model, Max, sum(y[c]*U[c] for c in C))
    
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

    for c in C
        # constraint for the variable Delta
        cons_y = @constraint(model, y[c]<=1)
        set_name(cons_y, "cons_delta_"*string(c))

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
            cons = @constraint(model, y[c] <= x[i, j]*ksi[i, j])
            set_name(cons, "cons_lambda_"*string(c)*"_"*string(i)*"_"*string(j))
        end 
    end
    
    for (v, C_v) in vertic_cycles
        # each node must be at most in one cycle
        cons = @constraint(model, sum(y[c] for c in C_v)<=1)
        set_name(cons, "cons_mu_"*string(v))
    end;

    return model
end
;
