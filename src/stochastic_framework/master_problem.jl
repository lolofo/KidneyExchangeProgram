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
    updateMasterProblem

Update the master problem to have the classical recourse modelisation

# Parameters
* `model` : a master problem type
* `duals` : the duals solutions associated with each scenarios (of shape |C| x S)
* `T_array` : the different scenarios (of shape |C| x S)
* `C` : our index cycles
* `U` : the utility of each cycle in our cycles (of shape |C|)
"""
function updateMasterProblem(model, duals, T_array, C, U)
    # TODO
    return Nothing
end
;





"""
    updateCVARMasterProblem

Update the master problem to have an risk adverse modelisation based on the CVAR criterion

# Parameters
* `model` : a master problem type
* `duals` : the duals solutions associated with each scenarios (of shape |C| x S)
* `T_array` : the different scenarios (of shape |C| x S)
* `C` : our index cycles
* `U` : the utility of each cycle in our cycles (of shape |C|)
"""
function updateCVARMasterProblem(model, duals, T_array, C, U)
    # TODO
    return Nothing
end
;