"""
This file will contain the algorithms we use to solve our different problems.
"""

using JuMP;
using GLPK;

include(join(["recourse_problem.jl"], Base.Filesystem.path_separator))
include(join(["master_problem.jl"], Base.Filesystem.path_separator))


"""
    LshapeClusterMethod

This method propose the Benders method for the clustering problem.
This method uses functions available in the files :

- recourse_problem.jl
- master_problem.jl

# Parameters
# `kep_graph` : the graph of the kidney exchange program
# `ClusterSize` : the size of the clusters
# `C` : the list of the cycle index
# `cycles` : the list of the real cycles
# `U` : the list of the utilities
# `vertic_cycles` : a Julia dictionnary. for each key the value is a list corresponding to the cycles involving the key
# `ksi` : the tensor of the scenarios
# `itmax` : the number of maximum iteration
# `verbose` : if true, the main steps will be printed on the standard output.

# Returns
This method returns a dictionnary with the following keys :
# `master_problem` : the master problem ready to be optimize.
# `first_level_var` : the x value.
# `objective_value` : the objective value.
"""
function LshapeClusterMethod(kep_graph, ClusterSize, C, cycles, U, vertic_cycles, ksi, itmax=100000, tol=1e-4, verbose = true)

    verbose && println("Start of the L-shape method for the cluster problem");

    # init of some variables
    master_val = 0

    # definition of the master problem
    res_master_problem = masterClusterProblem(kep_graph, ClusterSize, C, cycles, U, vertic_cycles)
    master_problem = res_master_problem["model"]

    # solve and store the solution
    optimize!(master_problem)
    sol = value.(master_problem[:x])
    
    # add the theta variables to our problem
    addThetaCluster(master_problem, size(ksi)[3], C, U);
    it = 0

    # init the different recourse problems
    recourse_pbs = []
    for k in 1:1:size(ksi)[3]
        model = recourseClusterProblem(sol, ksi[:, :, k], C, vertic_cycles, U, cycles)
        append!(recourse_pbs, [model])
    end

    for k in 1:1:size(ksi)[3]
        curr_recourse = recourse_pbs[k] # get the problem
        modifyRecourseClusterProblem(curr_recourse, sol, C, ksi[:, :, k]) # update with the scenar
        res_recourse = solveRecourseClusterProblem(curr_recourse, sol, C, vertic_cycles, cycles)
        dual = res_recourse["dual"] # get the dual solution
        updateCluster(kep_graph, master_problem, C, ksi[:, :, k], dual, k); # update the problem with the dual solution.
    end

    condition = true

    while condition

        # update the master problem and its solution
        optimize!(master_problem);
        sol = value.(master_problem[:x])
        master_val = objective_value(master_problem) # objective value of the master problem

        # solve the different recourse problems
        recourse_val = 0

        for k in 1:1:size(ksi)[3]
            curr_recourse = recourse_pbs[k] # get the problem
            modifyRecourseClusterProblem(curr_recourse, sol, C, ksi[:, :, k]) # update with the scenar
            res_recourse = solveRecourseClusterProblem(curr_recourse, sol, C, vertic_cycles, cycles)
            dual = res_recourse["dual"] # get the dual solution
            updateCluster(kep_graph, master_problem, C, ksi[:, :, k], dual, k); # update the problem with the dual solution.
            recourse_val += objective_value(curr_recourse)
        end
        
        recourse_val /= size(ksi)[3]

        ##################################################
        # the updates and the print along the iterations #
        ##################################################

        it += 1

        if it % 10 == 0
            verbose && print("Iteration ("*string(it)*") >> ") ;
            verbose && print("objective value : "*string(master_val))
            verbose && println()
        end
        
        # the stopping criterions
        if it >= itmax
            verbose && println("The maximum number of iteration is reached : stop");
            condition = false
        end

        # stopping criterion of
        if master_val <= (recourse_val + tol)
            verbose && println("The stopping criterion is reached");
            verbose && println("The objective value : "*string(master_val))
            condition = false
        end
    end
    return(Dict(
        "master_problem" => master_problem, 
        "first_level_var"=>sol,
        "objective_value" => master_val
    ))
end
;