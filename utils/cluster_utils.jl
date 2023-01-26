"""
    rank_index_cycle

rank cycles with respect to the probability and utility

# Parameters
*`data` : result from extractCycleInformation
"""
function rank_index_cycle(data)
    return reverse(sortperm(data["P"].*data["U"]))
end;