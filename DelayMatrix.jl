#using Graphs
function DelayMatrix(filename_Line)
    # filename_Line = "Line.csv"
    # branchmat = readtable(filename_Line)
    branchmat = readcsv(filename_Line, header=true)[1]
    fbus = Array{Int64}(1,size(branchmat,1))
    tbus = Array{Int64}(1,size(branchmat,1))
    for i = 1:size(branchmat,1)
        fbus[i] = branchmat[i,2]
        tbus[i] = branchmat[i,1]
    end
    N_bus = maximum(fbus)
    g = simple_inclist(N_bus)
    g = simple_inclist(N_bus,is_directed=false)
    dists = zeros(N_bus)
    for i = 1:size(branchmat,1)
        add_edge!(g,tbus[i],fbus[i])
        dists[i] = 1
    end
    delay_connected = zeros(N_bus,N_bus)
    for i = 1:N_bus
        temp = dijkstra_shortest_paths(g, dists, i)
        delay_connected[:,i] = temp.dists
    end

    delay_not_connected = zeros(N_bus,N_bus)
    for i = 1:N_bus
        for j = 1:N_bus
            delay_not_connected[i,j] = abs(i-j)
        end
    end

    delay = Array{Int64}(N_bus,N_bus)
    delay[:] = min.(delay_connected[:],delay_not_connected[:])

    return delay
end
