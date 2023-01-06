using JuMP;
using GLPK;

# the cluster problem
#####################

"""
recourseClusterProblem

This function will create the recourse problem with cycle formulation for multiple scenarios

# Parameters
* `x` : solution au MP
* `ksi` : random variable giving the success/fail of the cross test (one realisation of the random variable)
* `C` : list of cycle index
* `vertic_cycles` : Dictionary with {"graph vertices" : "lists of cycle which implie the keys"}
* `U` : list of utility of each cycle
* `cycle` : an exaustive list of the cycle

# Returns :
This functions will return
* `model` : the model of the recourseProblem
"""
function recourseClusterProblem(x, ksi, C, vertic_cycles, U, cycles)
    
    model = Model(GLPK.Optimizer)

    # rescouse variables (with relaxation)
    @variable(model, y[c in C] >= 0)
    @objective(model, Max, sum(y[c]*U[c] for c in C))

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
    end
    return(model)
end


"""
    solveRecourseClusterProblem

The previous function was made to create the recourse problem.
This new function will solve the recourse problem and get the dual of the problem

# Parameters
# `model` : the recourse problem
* `x` : solution au MP
* `ksi` : random variable giving the success/fail of the cross test (one realisation of the random variable)
* `C` : list of cycle index
* `vertic_cycles` : Dictionary with {"graph vertices" : "lists of cycle which implie the keys"}
* `U` : list of utility of each cycle
* `cycle` : an exaustive list of the cycle

# Returns :
This functions will return
* the duals solutions : the dual solution will allow us to update the master problem in the Benders decomposition.
The dual solutions are decomposed into three groups lambda, mu and Delta which corresponds to the different constraints of the recourse
problem.
"""
function solveRecourseClusterProblem(model, x, C, vertic_cycles, cycles)
    # optimisation 
    # handle exception because julia error aren't good.
    try
        optimize!(model)
    catch e
        println("optimization error : problem during the optimization part")
        println(e.msg)
    end
    ### Now lets recover all the duals solutions ###

    # the different dual solutions.
    dual_lambda = zeros(size(x)[1], size(x)[2], length(C))
    dual_mu = zeros(size(x)[1])
    dual_delta = zeros(length(C))
    
    if has_duals(model)
        for c in C
            # the dual for maximization is define as min -objectif so we have to take the opposite
            dual_delta[c] = -dual(constraint_by_name(model, "cons_delta_"*string(c)))
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
                # the dual for maximization is define as min -objectif so we have to take the opposite
                dual_lambda[i,j,c] = -dual(constraint_by_name(model, "cons_lambda_"*string(c)*"_"*string(i)*"_"*string(j)))
            end 
        end
        
        for (v, C_v) in vertic_cycles
            # the dual for maximization is define as min -objectif so we have to take the opposite
            dual_mu[v] = -dual(constraint_by_name(model, "cons_mu_"*string(v)))
        end
    else
        print("No dual : error stop the L shape methode")
        return(Nothing) # return nothing to raise exception after.
    end

    return Dict(
        "dual" => Dict("dual_lambda" => dual_lambda, "dual_mu" => dual_mu, "dual_delta" => dual_delta)
    )
end
;


"""
    modifyRecourseClusterProblem

This function will modify the recourseClusterProblem for the L-shape method.
The objective is, at each iteration of the L-shape method, to not create again another recourse problem.
There is only one constraint that need to be modified for our problem.

# Parameters
# `model` : the model of the recourse problem
# `x` : the first level value variable
# `C` : the cycle's indexes
# `ksi` : the scenario we are dealing with
"""
function modifyRecourseClusterProblem(model, x, C, ksi)

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

            # delete the constraint
            delete(model, constraint_by_name(model, "cons_lambda_"*string(c)*"_"*string(i)*"_"*string(j)))

            # add the new constraint
            cons = @constraint(model, model[:y][c] <= x[i, j]*ksi[i, j])
            set_name(cons, "cons_lambda_"*string(c)*"_"*string(i)*"_"*string(j))
        end 
    end
end;