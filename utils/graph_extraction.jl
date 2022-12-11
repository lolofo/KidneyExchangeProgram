using Graphs;


##############################################################################################
##############################################################################################
##############################################################################################

### Extract information for the cycles ###

"""
    getEdgeQuantiy

This function allow us to calculate the Edge probability and 
weight of a given cycle

# Parameters

*`g` : the kep graph
*`cycle` : a cycle in our kep_graph 
"""
function getEdgeQuantity(g, cycle)
    P = []
    W = []
    for v in [1:1:length(cycle);]
        s = cycle[v]
        d = Nothing
        if v < length(cycle)
            d = cycle[v+1]
        else
            d = cycle[1]
        end
        e = Edge((s, d))
        append!(P, get_prop(g, e, :failure))
        append!(W, get_prop(g, e, :weight))
    end
    return(Dict("P" => P, "W" => W))
end

"""
    getProbCycle

# Parameters
* `probs` : the probability of each edge in a cycle
"""
function getProbCycle(probs)
    buff_prob = 1
    for p in probs
        buff_prob *= (1 - p)
    end
    return(1 - buff_prob)
end


function getCycleUtility(weights, mode)
    # mode --> prod
    @assert mode ∈ ["prod", "sum", "max"] "the string mode must be in [prod, sum, max]"
    
    weight = Nothing

    # mode --> sum
    if mode == "prod"
        weight = 1
        for w in weights
            weight *= w
        end
    # mode --> sum
    elseif mode == "sum"
        weight = 0
        for w in weights
            weight += w
        end
    # mode --> max
    elseif mode == "max"
        weight = maximum(weights)
    end
    return(weight)
end;

""" 
    extractCycleInformation

This function allow us to extract the relevant information from the kep_graph.
The relevant information here concern the cycle formulation of the kep_graph.

# Parameters
* `g` : the kep_graph
* `K` (int): the length of the Cycles
* `mode`: the method to use to compute the calculus of the utilities

# Return 
This function returns a Julia dictionnary with the following keys :
* `Cycles_index` : a list of integer each element of the list corresponds to the index of a cycle
* `vertic_cycles` : a dictionnary with vertices as keys and a list of cycles which involve the key as value
* `Cycles` : the exhaustive enumeration of the cycles
* `P` : for each cycle, the probability of failure. To get the success do 1 - ...
* `U` : the utility of each cycle
"""
function extractCycleInformation(g, K, mode)

    enum_cycles = simplecycles_limited_length(g, K, 10^6)
    if length(enum_cycles)==0
        # in this case we return nothing
        println("Infeasible Problem")
        return(Nothing)
    end

    C = [1:1:length(enum_cycles);]

    # probabilities for each cycle to fail
    P = []
    U = []

    for c in C
        cycle = enum_cycles[c]
        res_edge = getEdgeQuantity(g, cycle)

        append!(P, getProbCycle(res_edge["P"]))
        append!(U, getCycleUtility(res_edge["W"], mode))
    end

    
    vertic_cycles = Dict()
    for v in vertices(g)
        C_v = []
        for c in C
            if v in enum_cycles[c]
                append!(C_v, c)
            end
        end
        merge!(vertic_cycles, Dict(v => C_v))
    end

    res = Dict("Cycles_index" => C, 
    "vertic_cycles" => vertic_cycles, 
    "Cycles" => enum_cycles, 
    "P" => P,
    "U" => U)

    return(res)
end
;


##############################################################################################
##############################################################################################
##############################################################################################

### Extract information for the cluster problem ###





