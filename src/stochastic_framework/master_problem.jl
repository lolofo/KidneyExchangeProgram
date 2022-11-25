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

    @objective(model, Min, sum( - x[c]U[c] for c in C))

    for (v, C_v) in vertic_cycles
        @constraint(model, sum(x[c] for c in C_v)<=1)
    end
    
    return(Dict("model" => model))

end
;

"""
TODO : faire une fonction qui ajoute au problème maitre 
       les contraintes en L (Benders)
"""