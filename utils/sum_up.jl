using DataFrames;
include(join(["solution_extraction.jl"], Base.Filesystem.path_separator))


"""
sum_up

This function is used to summarize instance dans display features such as quantity of nodes, edges before and after preprocessing. 
     
"""
function sum_up_instance(K)
    j = 1
    list_dir = readdir(join(["_cache", "data"], Base.Filesystem.path_separator));
    myregex = r"MD-00001-0+(\d+).dat"
    list_name = Int[]
    list_n_nodes = Int[]
    list_n_edges = Int[]
    list_n_circles = Int[]
    list_n_nodes_preprocess = Int[]
    list_n_edges_preprocess = Int[]
    
    for i in 1:1:(length(list_dir)/2)
        try
            kep_graph = read_kep_file("./_cache/data/"*list_dir[j+1],"./_cache/data/"*list_dir[j]);
            append!(list_name, parse(Int64, match(myregex, list_dir[j])[1]))
            append!(list_n_edges, ne(kep_graph))
            append!(list_n_nodes, nv(kep_graph))
            append!(list_n_circles, length(simplecycles_limited_length(kep_graph, K, 10^6)))
            # preprocessing 
            kep_graph, temp = removeUselessNodes(kep_graph, K)
            append!(list_n_edges_preprocess, ne(kep_graph))
            append!(list_n_nodes_preprocess, nv(kep_graph))
            j += 2
        catch
            println("error in reading")
            j = j+1
        end
    end

    return DataFrame(name = list_name, 
    n_edge=list_n_edges, 
    n_node=list_n_nodes, 
    n_circle=list_n_circles, 
    n_edge_preprocess=list_n_edges_preprocess, 
    n_node_preprocess=list_n_nodes_preprocess)
end;

"""
sum_up

This function pick one instance compute solution and evaluate the results. Finally, all resultsare stored in a dataframe.
     
# Parameters

*`number_instance` : number of instance (ex 1 -> MD-00001-00000001)
*`df` : dataframe which store the results
*`ClusterSize` : Cluser size
*`nb_scenar` : number of scenarios for the L-shaped
*`nb_scenar_eval` : number of scenarios for the solution evalution
*`nb_cycles` : quantity of cycle considered
*`K` : cycle size
*`dist` : distribution to use
*`cvar` : if true we solve problem with Cvar otherwise with expected value
*`risk_level` : risk level for Cvar
"""
function sum_up(number_instance, df, ClusterSize, nb_scenar, nb_scenar_eval, nb_cycles, K=2, dist="Constant", cvar = false, risk_level = 0)
    # read failure rate
    data = read_and_preprocess(number_instance, K, dist, nb_cycles)
    
    kep_graph = data["kep_graph"]
    
    if nv(kep_graph) == 0
        if cvar
            push!(df, [false 0 0 0 -1.0 "TODO" "no feas" "no feas" "None" "no feas" "no feas" "no feas" "no feas"])
        else
            push!(df, [false 0 0 0 -1.0 -1.0 "no feas" "no feas" "no feas" "no feas" "no feas" "no feas" "no feas"])
        end
    else
        C = data["Cycles_index"]
        cycles = data["Cycles"]
        U = data["U"]
        vertic_cycles = data["vertic_cycles"]

        # create scenarios
        ksi = getScenarioClusterK(kep_graph, nb_scenar)

        # solve problem
        res_mean = masterClusterProblem(kep_graph, ClusterSize, C, cycles, U, vertic_cycles)
        optimize!(res_mean["model"])

        res_lshape = LshapeClusterMethod(kep_graph, ClusterSize, C, cycles, U, vertic_cycles, ksi, 50, 1e-3, false, cvar, risk_level)
        
        obj_unroll = nothing # initialization
        if cvar
            # TODO implement unroll CVar problem
            obj_unroll = "TODO"
        else
            res_unroll = unrollClusterProblem(kep_graph, ClusterSize, C, cycles, U, ksi, vertic_cycles)
            optimize!(res_unroll["model"])
            obj_unroll = objective_value(res_unroll["model"])
        end
        
        if sum(value.(res_lshape["first_level_var"]))<0.1

            z_sp = "empty"
            z_ev = "empty"
            z_ws = "empty"
            VSS = "empty"
            EVPI = "empty"
            nb_pers_cluster = "empty"
            nb_cluster = "empty"
            nb_pers_transp = "empty"
            nb_cycle_transp = "empty"
        else

            res_z_sp = evaluateSolution_ls(kep_graph, nb_scenar_eval, res_lshape["first_level_var"], C, vertic_cycles, U, cycles)
            z_sp = res_z_sp["z_sp"]
            if cvar
                VSS = "None"
            else
                z_ev = evaluateSolution_ls(kep_graph, nb_scenar_eval, value.(res_mean["model"][:x]), C, vertic_cycles, U, cycles)["z_sp"]
                VSS = z_sp - z_ev
            end
            
            z_ws = evaluateSolution_ws(kep_graph, nb_scenar_eval, C, vertic_cycles, U, cycles, ClusterSize)
            # solution evaluation 
            
            EVPI = z_ws - z_sp

            nb_pers_cluster = res_z_sp["nb_pers_cluster"]
            nb_cluster = res_z_sp["nb_cluster"]
            nb_pers_transp = res_z_sp["nb_pers_transp"]
            nb_cycle_transp = res_z_sp["nb_cycle_transp"]
        end

            

        push!(df, 
            [res_lshape["optimal"]                # status of the problem
            res_lshape["nb_iterations"]           # nb benders iterations
            res_lshape["nb_added_constraints"]    # nb added benders contraints
            length(C)                             # number of cycles in the graph
            res_lshape["objective_value"]         # objective value at the end of the L-shape
            obj_unroll                            # unroll objective value
            z_sp                                  # stochastic objective value
            EVPI
            VSS
            nb_pers_cluster
            nb_cluster
            nb_pers_transp
            nb_cycle_transp])
    end
end;