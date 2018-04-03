using JuMP, Gurobi, Graphs, JLD, DataFrames
include("DataImport.jl")
include("TypeInfo.jl")
include("contin_data.jl")
include("load24_data.jl")
include("BFsub.jl")
include("DelayMatrix.jl")
#### Optimization Parameters Control Area ####
gap = 1e-2
delayMode = 1
##############################################
########## Simulation Case Setting ##########
N_line = 14
N_time = 24
N_storage = 3
N_bus = 15
N_scen_normal = 1
N_scen_contin = 1
N_scen = N_scen_normal + N_scen_contin
N_event = 1
##############################################
filename_Node = "Node.csv"
filename_Generator = "Generator.csv"
filename_Line = "Line.csv"
filename_Storage = "Storage1MWhx3.csv"
filename_Station = "Station.csv"
filename_SMP = "DataMiner-Export_2018-01-10-121427.csv"
filename_Load24 = "load24data-PJM.csv"

load24 = load24_data(filename_Load24)
# N_days = size(load24,1)



contin = contin_data(N_line,N_time,N_scen_normal,N_scen_contin,N_event)
scen_data = ScenarioData(contin,load24[1:N_scen_normal,:])

status,obj_value,simulation_time,ESlocation,SoC,u,z,v,a,pg,qg,IC,OC,EC,OCinter,OCgen,OCloadshed,OCstorage,x,sigma_d,sigma_l = BFsub(scen_data,gap,delayMode,filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP);

@save "result.jld"
