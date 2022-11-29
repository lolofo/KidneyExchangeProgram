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




