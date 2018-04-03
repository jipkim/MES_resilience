using DataFrames  # readtable function
# using Requests
type Bus
   nodeID::Int
   root::Int
   Pd::Float64
   Qd::Float64
   Vmax::Float64
   Vmin::Float64
   B::Float64
   R::Float64
   X::Float64
   children::Vector{Int}
   ancestor::Vector{Int}
   genids::Int
   function Bus(nodeID, root, Pd, Qd, B, R, X, Vmax, Vmin)
      b = new(nodeID, root, Pd, Qd)
      b.Vmax = Vmax
      b.Vmin = Vmin
      b.R = R
      b.X = X
      b.B = B
      b.children = Int[]
      b.ancestor = Int[]
      # b.genids = 0
      return b
   end
end
#################################################################
type Generator
   genID::Int
   busidx::Int
   Pgmax::Float64
   Pgmin::Float64
   Qgmax::Float64
   Qgmin::Float64
   cost::Float64
   function Generator( busidx, Pgmax, Pgmin, Qgmax, Qgmin, cost)
      g = new(busidx)
      g.busidx = busidx
      g.cost = cost
      g.Pgmax = Pgmax
      g.Pgmin = Pgmin
      g.Qgmax = Qgmax
      g.Qgmin = Qgmin
      return g
   end
end

##################################################################
type Line
   arcID::Int
   tail::Int # the "to" node
   head::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   u::Float64 # the capacity of the line
   function Line(arcID, tail, head, r, x, u)
      line = new(arcID, tail, head, r, x)
      line.u = u
      return line
   end
end
#########################################
## storage data
type Storage
   ID::Int
   busidx::Int
   Pkmax::Float64
   Ekmax::Float64
   Aleph::Float64
   Co::Float64
   Cp::Float64
   Ce::Float64
   function Storage(ID, busidx, Pkmax, Ekmax, Aleph, Co, Cp, Ce)
      s = new(ID)
      s.ID = ID
      s.busidx = busidx
      s.Pkmax = Pkmax
      s.Ekmax = Ekmax
      s.Aleph = Aleph
      s.Co = Co
      s.Cp = Cp
      s.Ce = Ce
      return s
   end
end
#########################################
## station data
type Station
   ID::Int
   busidx::Int
   Cf::Float64
   function Station(ID, busidx, Cf)
      st = new(ID)
      st.ID = ID
      st.busidx = busidx
      st.Cf = Cf
      return st
   end
end
######################################################################################################


function DataImport(filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP)

   ######################################################################################################
   # Bus/Node Data
   # busmat = readtable(filename_Node)
   busmat = readcsv(filename_Node, header=true)[1]

   buses = Bus[]
   busIDmap = Dict()
   for i in 1:size(busmat,1)

      nodeID = i
      busIDmap[busmat[i,1]] = i

      if i==1
         root = busIDmap[busmat[i,1]]
      else
         root = 0
      end

      Pd = busmat[i,2]
      Qd = busmat[i,3]
      R = busmat[i,6]
      X = busmat[i,7]

      Vmax = busmat[i,4]
      Vmin = busmat[i,5]
      B = busmat[i,8]

      b = Bus(nodeID, root, Pd, Qd, B, R, X, Vmax, Vmin)
      push!(buses, b)
   end

   #######################################################################################################
   ## generator data
   generatorlist = Int[]
   generators = Generator[]
   # genmat = readtable(filename_Generator)
   genmat = readcsv(filename_Generator, header=true)[1]

   for i in 1:size(genmat,1)
      busidx = genmat[i,1]
      Pgmax = genmat[i,2]
      Pgmin = genmat[i,3]
      Qgmax = genmat[i,4]
      Qgmin = genmat[i,5]
      cost = genmat[i,6]

      g = Generator(busidx, Pgmax, Pgmin, Qgmax, Qgmin, cost)
      push!(generators, g)
      # setg(buses[busidx], i)
   end
   #for g in 1:length(generators)
   #   generators[g].cost = genmat[g,4]
   #end
   #################################################################
   ## branch data
   # branchmat = readtable(filename_Line)
   branchmat = readcsv(filename_Line, header=true)[1]


   lines = Line[]
   for i in 1:size(branchmat,1)
      fbus = busIDmap[branchmat[i,2]]
      tbus = busIDmap[branchmat[i,1]]
      abus = busIDmap[branchmat[i,1]]
      x = branchmat[i,4]
      r = branchmat[i,3]
      u = branchmat[i,6]  # flow limit
      push!(buses[tbus].children, fbus)#children
      push!(buses[fbus].ancestor, abus)#ancestor
      l = Line(i, tbus, fbus, r, x, u)
      push!(lines,l)
   end
   #########################################
   ## storage data
   # storagemat = readtable(filename_Storage)
   storagemat = readcsv(filename_Storage, header=true)[1]

   storages = Storage[]
   for i in 1:size(storagemat,1)
      storageID = i
      busidx = 1
      Pkmax = storagemat[i,2]
      Ekmax = storagemat[i,3]
      Aleph = storagemat[i,4]
      Co = storagemat[i,5]
      Cp = storagemat[i,6]
      Ce = storagemat[i,7]

      temp = Storage(storageID, busidx, Pkmax, Ekmax, Aleph, Co, Cp, Ce)
      push!(storages,temp)
   end

   #########################################
   ## charging/discharing station data
   # stationmat = readtable(filename_Station)
   stationmat = readcsv(filename_Station, header=true)[1]
   stations = Station[]
   for i in 1:size(stationmat,1)
      stationID = i
      busidx = stationmat[i,2]
      Cf = stationmat[i,3]

      temp = Station(stationID, busidx, Cf)
      push!(stations,temp)
   end

   #########################################
   SMP_raw = readcsv(filename_SMP, header=true)[1]
   SMP = SMP_raw[3:end]

   return buses, generators,lines, storages, stations, SMP
end
