"""
    getScenario

This function return a scenario for the recourse problem

# Parameters
* `failureRate` : list of the probability that each cycle fails
"""
function getScenario(failureRate)
    T = [] 
    for fail in failureRate
        d = Binomial(1, fail)
        append!(T, rand(d, 1))
    end
    return(T)
end



"""
TODO : 
"""