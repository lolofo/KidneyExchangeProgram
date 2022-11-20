using JuMP;
using GLPK;

"""
    function to provide the recourse problem
"""
function recourseProblem(x, T, C, U)

    model = Model(GLPK.Optimizer)

    # relaxation on the recourse variable
    @variable(model, 1 >= y[i in C] >= 0)

    # 
    @objective(model, Min, sum(y[c]U[c] for c in C))

    @constraints(model, 
    begin
        recours_cons[c in C], y[c] = T[c]x[c]
    end
    )

    if has_duals(model)
        optimize!(model) # solve
        dual_solution = [dual(recours_cons[c]) for c in C]

    else 
        dual_solution = Nothing
        println("Error : no dual")

    return Dict(
        "model" => model,
        "dual" => dual_solution
    )
end


