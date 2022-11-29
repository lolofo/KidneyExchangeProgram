"""
    read_kep_file

Contruct a kep_graph from a `.wmd` and a `.dat` files from PrefLib.

# Parameters
* `wmd_file::String` : path of the `.wmd` file.
* `dat_file::String` : path of the `.dat` file.
"""
function read_kep_file(wmd_file::String, dat_file::String)
    wmd_file_name = split(split(wmd_file, '/')[end], '.')[1]
    dat_file_name = split(split(dat_file, '/')[end], '.')[1]

    wmd_file_name == dat_file_name || throw(ArgumentError(".wmd and .dat files do not correspond to the same dataset."))
    isfile(wmd_file) || throw(ArgumentError("$(wmd_file): file not found."))
    isfile(dat_file) || throw(ArgumentError(".dat file not found."))

    # Extract meta information from the .dat file
    file = readdlm(dat_file, '\n')
    nb_vertices = length(file)-1
    kep_graph = MetaDiGraph(nb_vertices, 0)
    for line in file[2:end]
        splitted_line = split(line, ',')
        pair = parse(Int, splitted_line[1])
        set_prop!(kep_graph, pair, :pra, parse(Float64, splitted_line[5]))
    end
    # Extract the graph structure from the .wmd file using a MetaGraph
    # get the number of vertices in first line
    wmd_io = open(wmd_file, "r")

    # skip next nb_vertices lines, which are redundant with the data contained in the .dat file
    for i in 1:nb_vertices+1
        readline(wmd_io)
    end
    # read the set of edges
    while !eof(wmd_io)
        splitted_line = split(readline(wmd_io), ',')
        # /!\ Pairs are numbered from 0 in the second part of the file
        source = parse(Int, splitted_line[1]) + 1
        destination = parse(Int, splitted_line[2]) + 1
        weight = parse(Int, splitted_line[3])

        # do not add an edge that has a zero weight
        if weight > 0
            add_edge!(kep_graph, source, destination, :weight, weight)
        end
    end
    return kep_graph
end
;

################################################################################################################
################################################################################################################
################################################################################################################

DISTRIBUTIONS = ["Constant","Binomial","BinomialUNOS","BinomialAPD","NoFailure"]

"""
    get_failure_rates

Generate failure rates on each edge, and add its value as a property to the edge of the kep_graph. 

# Parameters
* `kep_graph::MetaDiGraph` : graph describing the pairs and compatibilities
* `distribution::String` : type of distirbution of uncertainties; to be chosen in the DISTRIBUTIONS vector
"""
function get_failure_rates(kep_graph::MetaDiGraph, distribution::String)

    failure_rates = []

    for edge in edges(kep_graph)
        # Failure rates depend on the chosen distribution of uncertainties
        if distribution == "Constant"
            # constant failure rates equal to 70%
            set_prop!(kep_graph, edge, :failure, 0.21)
        elseif distribution == "Binomial"
            if rand() < 0.25
                # random failure rates equal to 10% on average for 25% edges
                set_prop!(kep_graph, edge, :failure, rand() * 0.2)
            else
                # random failure rates equal to 90% on average for 75% edges
                set_prop!(kep_graph, edge, :failure, 0.8 + rand() * 0.2)
            end
        elseif distribution == "BinomialUNOS"
            # %pra denotes the panel reactive antibody level
            # %pra of the patient < 0.8 : UNOS low sensitized patients
            if get_prop(kep_graph, edge.dst, :pra) < 0.8
                # failure rate equal to 10% if the patient is low sensitized
                set_prop!(kep_graph, edge, :failure, 0.1)
            else
                # failure rate equal to 90% otherwise 
                set_prop!(kep_graph, edge, :failure, 0.9)
            end
        elseif distribution  == "BinomialAPD"
            # %pra denotes the panel reactive antibody level
            # %pra of the patient < 0.75 : APD low sensitized patients
            if get_prop(kep_graph, edge.dst, :pra) < 0.75
                # failure rate equal to 28% if the patient is low sensitized
                set_prop!(kep_graph, edge, :failure, 0.28)
            else
                # failure rate equal to 58% otherwise 
                set_prop!(kep_graph, edge, :failure, 0.58)
            end
        elseif distribution == "NoFailure"
            # failure rates equal to 0
            set_prop!(kep_graph, edge, :failure, 0.)
        end
    end
end
;