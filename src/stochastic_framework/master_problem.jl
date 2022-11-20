function masterProblem(C, vertic_cycles ,U)

    model = Model(GLPK.Optimizer)

    @variable(model, x[i in C], Bin)

    @objective(model, Min, sum( - x[c]U[c] for c in C))

    for C_v in vertic_cycles
        @constraint(model, sum(x[c] for c in C_v)<=1)
    end

    return Dict("model" => model)

end
;

"""
TODO : faire une fonction qui ajoute au problème maitre 
       les contraintes en L (Benders)
"""