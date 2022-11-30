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
        append!(T, rand(d, 1))
    end
    return(T)
end


function getScenarioK(failureRate, K)
    card_cycle = length(failureRate)
    T = zeros(Int8, card_cycle, K)
    for i in 1:1:K
        T[:, i] = getScenario(failureRate)
    end
    return(T)
end


;