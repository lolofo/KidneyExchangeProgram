using DataFrames;
include(join(["cluster_utils.jl"], Base.Filesystem.path_separator))

"""

"""
function transform_vertic_cycle(vertic_cycles, C)
    v_c = Dict()
    for k in vertic_cycles
        key = k[1]
        v_c[key] = []
        for c in vertic_cycles[key]
            if c ∈ C 
                append!(v_c[key], c)
            end;
        end;
        if length(v_c[k[1]])==0
            pop!(v_c, k[1])
        end;
    end;
    return v_c
end


"""

"""
function sum_up(kep_graph, df, ClusterSize, nb_scenar, nb_scenar_eval, nb_cycles, K=2, dist="Constant")
    # read failure rate
    failure_rates = get_failure_rates(kep_graph, dist);

    # read data
    data = extractCycleInformation(kep_graph, K, "sum");


    rank_index = rank_index_cycle(data) # les indexes des cycles dont on a besoin
    nb_cycles = min(nb_cycles, length(data["Cycles"]))

    C = data["Cycles_index"][rank_index][1:nb_cycles]
    cycles = data["Cycles"]
    U = data["U"]
    vertic_cycles = transform_vertic_cycle(data["vertic_cycles"], C)

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
        length(C)                             # number of cycles in the graph
        res_lshape["objective_value"]         # objective value at the end of the L-shape
        objective_value(res_unroll["model"])  # unroll objective value
        z_sp                                  # stochastic objective value
        EVPI
        VSS])

    return df
end;