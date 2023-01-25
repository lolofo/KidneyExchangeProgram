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

    # create scenarios
    ksi = getScenarioClusterK(kep_graph, nb_scenar)

    # solve problem
    res_mean = masterClusterProblem(kep_graph, ClusterSize, C, cycles, U, vertic_cycles)
    res_lshape = LshapeClusterMethod(kep_graph, ClusterSize, C, cycles, U, vertic_cycles, ksi, 50, 1e-3, false)
    res_unroll = unrollClusterProblem(kep_graph, ClusterSize, C, cycles, U, ksi, vertic_cycles)
    
    optimize!(res_mean["model"])
    optimize!(res_unroll["model"])

    z_sp = evaluateSolution_ls(kep_graph, nb_scenar_eval, res_lshape["first_level_var"], C, vertic_cycles, U, cycles)
    z_ev = evaluateSolution_ls(kep_graph, nb_scenar_eval, value.(res_mean["model"][:x]), C, vertic_cycles, U, cycles)
    z_ws = evaluateSolution_ws(kep_graph, nb_scenar_eval, C, vertic_cycles, U, cycles, ClusterSize)

    # solution evaluation 
    VSS = z_sp - z_ev
    EVPI = z_ws - z_sp

    push!(df, 
        [res_lshape["optimal"]                # status of the problem
        res_lshape["nb_iterations"]           # nb benders iterations
        res_lshape["nb_added_constraints"]    # nb added benders contraints
        length(cycles)                        # number of cycles in the graph
        res_lshape["objective_value"]         # objective value at the end of the L-shape
        objective_value(res_unroll["model"])  # unroll objective value
        z_sp                                  # stochastic objective value
        EVPI
        VSS])

    return df

end