"""
This file will contain the algorithms we use to solve our different problems.
"""


using JuMP;
using GLPK;

include(join(["recourse_problem.jl"], Base.Filesystem.path_separator))
include(join(["master_problem.jl"], Base.Filesystem.path_separator))

"""
    checkClusterConstraint

This function will check if the constraint provided by the dual is good enough to be added to the master

# Parameters
* `theta` : master variable theta
* `dual` : dictionnary of the dual variables
* `ksi` : the scenario
* `tol` : the tolerance to reach the optimality

# Returns
This function returns a boolean.
* `true` : we should update the master problem and add the constraint which corresponds to the dual
* `false` : the master problem should not be updated
"""
function checkClusterConstraint(theta, dual, x, ksi, tol)
    
    # our dual variables
    μ = dual["dual_mu"]
    λ = dual["dual_lambda"]
    δ = dual["dual_delta"]

    # the dual function
    s1 = sum(λ .* ksi .* x)
    s2 = sum(μ)
    s3 = sum(δ)
    s = s1 + s2 + s3 + tol

    # should we add the constraint
    # if theta < s ⇒ the constraint is already satisfied (useless)
    res = theta < s ? false : true

    return res
end;


"""
    LshapeClusterMethod

This method propose the Benders method for the clustering problem.
This method uses functions available in the files :

- recourse_problem.jl
- master_problem.jl

# Parameters
* `kep_graph` : the graph of the kidney exchange program
* `ClusterSize` : the size of the clusters
* `C` : the list of the cycle index
* `cycles` : the list of the real cycles
* `U` : the list of the utilities
* `vertic_cycles` : a Julia dictionnary. for each key the value is a list corresponding to the cycles involving the key
* `ksi` : the tensor of the scenarios
* `itmax` : the number of maximum iteration
* `verbose` : if true, the main steps will be printed on the standard output.
* `cvar`: 

# Returns
This method returns a dictionnary with the following keys :
* `first_level_var` : the value of x i.e. the first level solution.
* `objective_value` : the objective value.
* `nb_added_constraints` : the number of constraints we added through the iterations of the algorithm
* `optimal` : a boolean value : true : the problem is solved
* `nb_iterations` : the number of iterations the algorithm did
"""
function LshapeClusterMethod(
    kep_graph, 
    ClusterSize, 
    C, 
    cycles, 
    U, 
    vertic_cycles, 
    ksi,
    itmax=100000, tol=1e-4, verbose = true,
    cvar = false,
    risk_level = 0)

    verbose && println("Start of the L-shape method for the cluster problem");

    # __init__
    ##########
    nb_added_constraints = 0  # number of added constraint to the master problem
    it = 0                    # number of iterations
    master_val = 0            # the master problem objective value

    res_master_problem = masterClusterProblem(kep_graph, ClusterSize, C, cycles, U, vertic_cycles)
    master_problem = res_master_problem["model"]

    optimize!(master_problem)
    sol = value.(master_problem[:x])

    # update the master problem for our specific task
    if cvar
        addCVaRVariables(master_problem, size(ksi)[3], risk_level)
    else
        addThetaCluster(master_problem, size(ksi)[3], C, U);
    end;
    
    recourse_pbs = []
    for k in 1:1:size(ksi)[3]
        model = recourseClusterProblem(sol, ksi[:, :, k], C, vertic_cycles, U, cycles)
        append!(recourse_pbs, [model])
    end;

    for k in 1:1:size(ksi)[3]
        curr_recourse = recourse_pbs[k]
        modifyRecourseClusterProblem(curr_recourse, sol, C, cycles, ksi[:, :, k])
        res_recourse = solveRecourseClusterProblem(curr_recourse, sol, C, vertic_cycles, cycles)
        dual = res_recourse["dual"]
        updateCluster(kep_graph, master_problem, C, ksi[:, :, k], dual, k)
    end;

    # optimization loop 
    ####################

    optimal = false

    while !optimal
        it += 1

        # update the master problem and its solution
        optimize!(master_problem);
        sol = value.(master_problem[:x])                 # x value
        θ = value.(master_problem[:theta])               # θ values
        master_val = objective_value(master_problem)     # objective value

        cpt = 0

        # ∀ the scenarios
        for k in 1:1:size(ksi)[3]
            curr_recourse = recourse_pbs[k]
            modifyRecourseClusterProblem(curr_recourse, sol, C, cycles, ksi[:, :, k])
            res_recourse = solveRecourseClusterProblem(curr_recourse, sol, C, vertic_cycles, cycles)
            dual = res_recourse["dual"]
            # check if the dual got us a good constraint
            buff = checkClusterConstraint(θ[k], dual, sol, ksi[:, :, k], tol)
            if buff
                updateCluster(kep_graph, master_problem, C, ksi[:, :, k], dual, k)
                cpt += 1
                nb_added_constraints += 1
            end;
        end;

        # updates and log
        #################
        
        # logs every 10 iterations
        if it % 10 == 0
            verbose && print("Iteration ("*string(it)*") >> ") ;
            verbose && print("objective value : "*string(master_val))
            verbose && println()
        end;

        # no constraint added ⇒ the optimality is reached
        if cpt == 0
            optimal = true
        end;

        # maximum number of ietration reached
        if it >= itmax
            verbose && println("The maximum number of iteration is reached : stop");
            optimal = true
        end;

        # the optimality is reached
        if optimal
            verbose && println("The stopping criterion is reached");
            verbose && println("The objective value : "*string(master_val))
        end;
    end;

    # at last we check is the problem is at an optimal status
    ##########################################################

    optimal = it < itmax ? true : false
    !optimal && println("max it reached ⇒ increase the maximum number of iterations")

    return(Dict(
        "first_level_var"=>sol,
        "objective_value" => master_val,
        "nb_added_constraints" => nb_added_constraints,
        "optimal" => optimal,
        "nb_iterations" => it
    ));
end
;