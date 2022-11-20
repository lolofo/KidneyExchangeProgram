"""
Here we will define the first model for our problem
"""

using Graphs;
using JuMP;
using GLPK;




function classical_cycle_model(g, utilities, n)
    enum_cycles = simplecycles_limited_length(g, n, 10^6)
    if length(enum_cycles)==0
        # in this case we return nothing
        println("Infeasible Problem")
        return(Nothing)
    end

    # then we can start to create the model
    
    model = Model(GLPK.Optimizer)
    C = [1:1:length(enum_cycles);]

    # define the utility of each cycle
    U = []
    for c in C
        append!(U, sum(utilities[v] for v in enum_cycles[c]))
    end


    @variable(model, y[i in C], Bin) # one variable for each cycle
    @objective(model, Max, sum(y[c]U[c] for c in C))
    for v in vertices(g)
        C_v = []
        for c in C
            if v in enum_cycles[c]
                append!(C_v, c)
            end
        end
        @constraint(model, sum(y[c] for c in C_v)<=1)
    end
    return model, enum_cycles
end
;