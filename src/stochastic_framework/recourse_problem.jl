using JuMP;
using GLPK;

"""
    recourseProblem

This function will create the recourse problem 

# Parameters
* `x` : the value of the first level solution of shape |C|
* `T` : a given scenario
* `C` : the index of the cycles
* `U` : the utilities of each cycles
"""
function recourseProblem(x, T, C, U)

    model = Model(GLPK.Optimizer)

    # relaxation on the recourse variable
    @variable(model, 1 >= y[i in C] >= 0)

    @objective(model, Min, sum(y[c]U[c] for c in C))

    @constraints(model, 
    begin
        recours_cons[c in C], y[c] == T[c]x[c]
    end
    )

    # try to solve the recourse problem
    optimize!(model)

    if has_duals(model)
        dual_solution = [dual(recours_cons[c]) for c in C]
    else 
        dual_solution = Nothing
        println("Error : no dual")
    end

    return Dict(
        "model" => model,
        "dual" => dual_solution
    )
end;



"""
    solveRecourse

This function will solve the recourse problem for multiple scenarios

Some notations :
    |C| : the number of cycles
    S : the number of scenarios

# Parameters
* `x` : the first level solution (an array of binary variables)
* `T_array` : matrix of the different scenarios (of shape |C| x S)
* `C` : the index of each cycles in the kep graph (of shape |C|)
* `U` : the utility of each cycle in the kep_graph (of shape |C|)
"""
function solveRecourse(x, T_array, C, U)

    @assert length(x) == size(T_array)[1] "dimension problem between x and T_array"
    @assert length(x) == length(C) "dimension problem between x and C"
    @assert length(C) == length(U) "dimension problem between C and U"

    duals_solution = zeros(length(C), size(T_array))

    # iterate through the different scenarios
    for k in 1:1:size(T_array)[2]
        # select a scenario
        T = T_array[:, k]

        # solve the recourse problem corresponding to this scenario
        res = recourseProblem(x, T, C, U)

        # get the dual solution from this problem
        if res["dual"] != Nothing
            duals_solution[:,k] = res["dual"]
        end
    end

    # return the different dual values
    return Dict(
        "duals" => duals_solution
    )
end 
;

################################################################################################
################################################################################################
################################################################################################

"""
recourseClusterProblem

This function will solve the recourse problem with cycle formulation for multiple scenarios

# Parameters
* `x` : solution au MP
* `ksi` : random variable giving the success/fail of the cross test (one realisation of the random variable)
* `C` : list of cycle index
* `vertic_cycles` : Dictionary with {"graph vertices" : "lists of cycle which implie the keys"}
* `U` : list of utility of each cycle
* `cycle` : an exaustive list of the cycle

# Returns :
This functions will return
* the model : The model will allow us to access some other values later.
* the duals solutions : the dual solution will allow us to update the master problem in the Benders decomposition.
The dual solutions are decomposed into three groups lambda, mu and Delta which corresponds to the different constraints of the recourse
problem.
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
            dual_delta[c] = dual(constraint_by_name(model, "cons_y_"*str(c)))

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
                dual_lambda[i,j,c] = dual(constraint_by_name(model, "cons_lambda_"*string(c)*"_"*string(i)*"_"*string(j)))
            end 
        end
        
        for (v, C_v) in vertic_cycles
            dual_mu[v] = dual(constraint_by_name(model, "cons_mu_"*string(v)))
        end
    else
        print("No dual : error stop the L shape methode")
        return(Nothing) # return nothing to raise exception after.
    end

    return Dict(
        "model" => model,
        "dual" => Dict("dual_lambda" => dual_lambda, "dual_mu" => dual_mu, "dual_delta" => dual_delta)
    )
end
;

