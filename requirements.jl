using Pkg

dependencies = [
    "Random", 
    "MetaGraphs", 
    "SimpleWeightedGraphs", 
    "Graphs", 
    "JuMP", 
    "DelimitedFiles", 
    "Distributions",  
    "GraphPlot", 
    "HiGHS", 
    "DataFrames", 
    "Plots", 
    "BenchmarkTools",
    "Statistics"
]

Pkg.add(dependencies)