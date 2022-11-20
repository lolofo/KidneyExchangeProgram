using Graphs;
using JuMP;
using GLPK;

function likelihood_cycle_model(g, utilities, n)

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

    L = [] # log likelihood of each cycle

    for c in C
        append!(L, 0) # init the current likelihood
        for v in [1:1:length(enum_cycles[c]);]
            s = enum_cycles[c][v] # source of the edge
            d = Nothing
            if v < length(c)
                d = enum_cycles[c][v+1]
            else
                d = enum_cycles[c][1]
            end
            e = Edge.([(s, d)])
            # add 10^-16 for numerical stability
            L[c] = L[c] + log(1 - get_prop(g, e[1], :failure) + 10^(-16))
        end
    end


    @variable(model, y[i in C], Bin) # one variable for each cycle
    @objective(model, Max, sum(y[c]U[c] + y[c]L[c] for c in C))
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