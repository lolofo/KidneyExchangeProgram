using DataFrames;
include(join(["solution_extraction.jl"], Base.Filesystem.path_separator))



"""

"""
function sum_up(number_instance, df, ClusterSize, nb_scenar, nb_scenar_eval, nb_cycles, K=2, cvar=false, dist="Constant")
    # read failure rate
    data = read_and_preprocess(number_instance, K, dist, nb_cycles)
    
    kep_graph = data["kep_graph"]
    
    if nv(kep_graph) == 0
        push!(df, [false 0 0 0 -1.0 -1.0 -100 -100 -100 -100 -100 -100 -100])
    else
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

        res_z_sp = evaluateSolution_ls(kep_graph, nb_scenar_eval, res_lshape["first_level_var"], C, vertic_cycles, U, cycles)
        z_sp = res_z_sp["z_sp"]
        z_ev = evaluateSolution_ls(kep_graph, nb_scenar_eval, value.(res_mean["model"][:x]), C, vertic_cycles, U, cycles)["z_sp"]
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
            VSS
            res_z_sp["nb_pers_cluster"]
            res_z_sp["nb_cluster"]
            res_z_sp["nb_pers_transp"]
            res_z_sp["nb_cycle_transp"]])
    end
end;