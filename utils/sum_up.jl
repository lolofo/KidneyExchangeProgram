using DataFrames;



function sum_up(kep_graph, df, ClusterSize, nb_scenar, nb_scenar_eval)
    # read failure rate
    failure_rates = get_failure_rates(kep_graph, "Constant");

    # read data
    data = extractCycleInformation(kep_graph, 3, "sum");
    C = data["Cycles_index"] 
    cycles = data["Cycles"]
    U = data["U"]
    vertic_cycles = data["vertic_cycles"]

    println("data read")
    # create scenarios
    ksi = getScenarioClusterK(kep_graph, nb_scenar)
    println("create scenario")
    # solve problem
    res_mean = masterClusterProblem(kep_graph, ClusterSize, C, cycles, U, vertic_cycles)
    res_lshape = LshapeClusterMethod(kep_graph, ClusterSize, C, cycles, U, vertic_cycles, ksi, 50, 1e-3, false)
    res_unroll = unrollClusterProblem(kep_graph, ClusterSize, C, cycles, U, ksi, vertic_cycles)
    
    
    optimize!(res_mean["model"])
    optimize!(res_unroll["model"])

    
    # solution evaluation 
    EVPI = res_lshape["objective_value"] - evaluateSolution_ls(kep_graph, nb_scenar_eval, value.(res_mean["model"][:x]), C, vertic_cycles, U, cycles)
    println("create solve models")
    WS = evaluateSolution_ws(kep_graph, nb_scenar_eval, C, vertic_cycles, U, cycles, ClusterSize) - res_lshape["objective_value"]
    println(evaluateSolution_ws(kep_graph, nb_scenar_eval, C, vertic_cycles, U, cycles, ClusterSize))
    println(res_lshape["optimal"])

    push!(df, [res_lshape["optimal"] 
        res_lshape["nb_iterations"] 
        res_lshape["nb_added_constraints"] 
        res_lshape["objective_value"] 
        objective_value(res_unroll["model"]) 
        EVPI
        WS])

    return df

end