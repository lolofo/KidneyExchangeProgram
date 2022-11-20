using Graphs;

function extractInformation(g, K, compute_utility)

    enum_cycles = simplecycles_limited_length(g, K, 10^6)
    if length(enum_cycles)==0
        # in this case we return nothing
        println("Infeasible Problem")
        return(Nothing)

    C = [1:1:length(enum_cycles);]

    # probabilities for each cycle to fail
    P = []
    for c in C
        buff = 1
        for v in [1:1:length(enum_cycles[c]);]
            s = enum_cycles[c][v]
            d = Nothing
            if v < length(c)
                d = enum_cycles[c][v+1]
            else
                d = enum_cycles[c][1]
            end
            e = Edge.([(s, d)])
            buff *= (1 - get_prop(g, e[1], :failure))
        append!(P, 1 - buff)

    """
    TODO : gérer la manière dont on calcul les utilités pour
           cycles.
    """

    
    return Dict(
        "Cycles_index" => C,
        "Cycles" => enum_cycles,
        "U" => U,
        "P" => P
    )
end


