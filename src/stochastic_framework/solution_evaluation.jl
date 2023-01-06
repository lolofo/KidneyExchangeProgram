"""
In this file we propose a series of methods which will 

"""

include(join(["src", "stochastic_framework", "recourse_problem.jl"], Base.Filesystem.path_separator))
include(join(["src", "utils", "monte_carlo.jl"], Base.Filesystem.path_separator))


function evaluate_solution(nb_scenar, x, C, vertic_cycles, U, cycles)
    ksi = getScenarioClusterK(kep_graph, nb_scenar)
    j = 0
    obj = 0
    for i in 1:1:(nb_scenar)
        if j == 0
            model = recourseClusterProblem(x, ksi[:, :, i], C, vertic_cycles, U, cycles)["model"]
            obj += objective_value(model)
        else
            modifyRecourseClusterProblem(model, x, C, ksi[:, :, i])
            obj += objective_value(model)
        end
    end
    return obj / nb_scenar
end