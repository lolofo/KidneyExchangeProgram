"""
    getScenario

This function return a scenario for the recourse problem

# Parameters
* `failureRate` : list of the probability that each cycle fails (of shape |C|)
"""
function getScenario(failureRate)
    T = [] 
    for fail in failureRate
        d = Binomial(1, fail)
        append!(T, rand(d, 1)[1])
    end
    return T
end


"""
getScenarioK

This function return K scenarios for the recourse problem

# Parameters
* `failureRate` : list of the probability that each cycle fails (of shape |C|)
* `K` : scenario amount
"""

function getScenarioK(failureRate, K)
    card_cycle = length(failureRate)
    T = zeros(Int8, card_cycle, K)
    for i in 1:1:K
        T[:, i] = getScenario(failureRate)
    end
    return T
end



"""
getScenarioCluster

This function return a scenario for the recourse problem with cluster formulation

# Parameters
* `kep_graph` : instance graph
"""
function getScenarioCluster(kep_graph)
    n = nv(kep_graph)
    S = zeros(n, n)
    for i in 1:1:n
        for j in 1:1:n
            if has_edge(kep_graph, i, j)
                d = Binomial(1, 1-get_prop(kep_graph, Edge((i, j)), :failure))
                S[i, j] = rand(d, 1)[1]
            end 
        end
    end
    return S
end
;


"""
getScenarioK

This function return K scenarios for the recourse problem with cluster formulation

# Parameters
* `kep_graph` : instance graph
"""

function getScenarioClusterK(kep_graph, K)
    n = nv(kep_graph) 
    Sk = zeros(n, n, K) 
    
    for i in 1:1:K
        Sk[:, :, i] =  getScenarioCluster(kep_graph)
    end
    return Sk
end