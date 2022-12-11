using Random, MetaGraphs, SimpleWeightedGraphs, Graphs, JuMP, DelimitedFiles, Distributions, GLPK

""" 
    completRecourseProblem

# Parameters
* `C` : index of the cycles
* `vertic_cycles` : a dictionnary, at the key i of this list, give the a list of the index of the cycles which involve the node i
* `U` : the utilities of each cycle in the graph.
* `S_P` : probability of sucess of each cycle
"""
function completRecourseProblem(C, vertic_cycles , U, S_P)

    model = Model(GLPK.Optimizer)

    @variable(model, x[i in C], Bin)

    @objective(model, Max, sum(x[c]U[c]S_P[c] for c in C))

    for (v, C_v) in vertic_cycles
        @constraint(model, sum(x[c] for c in C_v)<=1)
    end
    
    return(Dict("model" => model))

end
;