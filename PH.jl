using JuMP, Gurobi, Graphs, JLD, DataFrames
include("DataImport.jl")
include("TypeInfo.jl")
include("contin_data.jl")
include("load24_data.jl")
include("PHsub_ini.jl")
include("PHsub_iter.jl")
include("DelayMatrix.jl")

#### Optimization Parameters Control Area ####
gap_ini = 1e-3
gap = 1e-3
rho_zero = 70.962 # 1MWh/300kW => 124.1835  #1MWh/150kW => 70.962
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
#### Progressive Hedging Penalty Setting ####
rho_x = rho_zero * ones(N_storage)
rho_z = (rho_zero / 1 )* ones(N_storage,N_bus)
rho = PenaltyParameter(rho_z,rho_x)
##############################################
filename_Node = "Node.csv"
filename_Generator = "Generator.csv"
filename_Line = "Line.csv"

filename_Storage = "Storage1MWhx3.csv"

filename_Station = "Station.csv"
filename_SMP = "DataMiner-Export_2018-01-10-121427.csv"
filename_Load24 = "load24data-PJM.csv"

load24 = load24_data(filename_Load24)

contin = contin_data(N_line,N_time,N_scen_normal,N_scen_contin,N_event)
scen_data = ScenarioData(contin,load24)

S_N = Int64[]
S_E = Int64[]
for s = 1:N_scen
    if sum(contin,1:2)[s] == 0
        push!(S_N,s)
    elseif sum(contin,1:2)[s] > 0
        push!(S_E,s)
    end
end

Prob = Array{Float64}(N_scen,1) # Scenario Probability
Prob[S_N] = 0.9 / length(S_N)
Prob[S_E] = 0.1 / length(S_E)


status_ini = Array{String}(N_scen,1)
obj_value_ini = Array{Float64}(N_scen,1)
simulation_time_ini = Array{Float64}(N_scen,1)
ESlocation_ini = Array{Int64}(N_storage,N_time,N_scen)
SoC_ini = Array{Float64}(N_storage, N_time, N_scen)
u_ini = Array{Float64}(N_storage, N_bus, N_time, N_scen)
z_ini = Array{Float64}(N_storage, N_bus, N_scen)
v_ini = Array{Float64}(N_bus, N_time, N_scen)
a_ini = Array{Float64}(N_bus, N_time, N_scen)
pg_ini = Array{Float64}(N_bus, N_time, N_scen)
qg_ini = Array{Float64}(N_bus, N_time, N_scen)
IC_ini = Array{Float64}(N_scen,1)
OC_ini = Array{Float64}(N_scen,1)
EC_ini = Array{Float64}(N_scen,1)
OCinter_ini = Array{Float64}(N_scen,1)
OCgen_ini = Array{Float64}(N_scen,1)
OCloadshed_ini = Array{Float64}(N_scen,1)
OCstorage_ini = Array{Float64}(N_scen,1)
x_ini = Array{Float64}(N_storage, N_scen)
sigma_d_ini = Array{Float64}(N_bus, N_time, N_scen)
sigma_l_ini = Array{Float64}(N_line, N_time, N_scen)


#

for s = 1:N_scen
    if s in S_N
        contin_temp = contin[:,:,s]
        scen_data_temp = ScenarioData(contin_temp,load24[s,:])
    elseif s in S_E
        contin_temp = contin[:,:,s]
        scen_data_temp = ScenarioData(contin_temp,mean(load24[S_N,:],1))
    end
    status_ini[s],obj_value_ini[s],simulation_time_ini[s],ESlocation_ini[:,:,s],SoC_ini[:,:,s],u_ini[:,:,:,s],z_ini[:,:,s],v_ini[:,:,s],a_ini[:,:,s],pg_ini[:,:,s],qg_ini[:,:,s],IC_ini[s],OC_ini[s],EC_ini[s],OCinter_ini[s],OCgen_ini[s],OCloadshed_ini[s],OCstorage_ini[s],x_ini[:,s],sigma_d_ini[:,:,s],sigma_l_ini[:,:,s]  = PHsub_ini(scen_data_temp, gap_ini, delayMode, filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP)
end

@save "PH_ini.jld"
#@load "PH_ini.jld"

#############
status = Array{String}(N_scen,1)
obj_value = Array{Float64}(N_scen,1)
simulation_time = Array{Float64}(N_scen,1)
ESlocation = Array{Int64}(N_storage,N_time,N_scen)
SoC = Array{Float64}(N_storage, N_time, N_scen)
u = Array{Float64}(N_storage, N_bus, N_time, N_scen)
z = Array{Float64}(N_storage, N_bus, N_scen)
v = Array{Float64}(N_bus, N_time, N_scen)
a = Array{Float64}(N_bus, N_time, N_scen)
pg = Array{Float64}(N_bus, N_time, N_scen)
qg = Array{Float64}(N_bus, N_time, N_scen)
IC = Array{Float64}(N_scen,1)
OC = Array{Float64}(N_scen,1)
EC = Array{Float64}(N_scen,1)
OCinter = Array{Float64}(N_scen,1)
OCgen = Array{Float64}(N_scen,1)
OCloadshed = Array{Float64}(N_scen,1)
OCstorage = Array{Float64}(N_scen,1)
x = Array{Float64}(N_storage, N_scen)
sigma_d = Array{Float64}(N_bus, N_time, N_scen)
sigma_l = Array{Float64}(N_line, N_time, N_scen)
#########################
N_iter = 40
#########################
status_iter = Array{String}(N_scen, N_iter)
obj_value_iter = Array{Float64}(N_scen, N_iter)
simulation_time_iter = Array{Float64}(N_scen, N_iter)
ESlocation_iter = Array{Int64}(N_storage,N_time,N_scen, N_iter)
SoC_iter = Array{Float64}(N_storage, N_time, N_scen, N_iter)
u_iter = Array{Float64}(N_storage, N_bus, N_time, N_scen, N_iter)
z_iter = Array{Float64}(N_storage, N_bus, N_scen, N_iter)
v_iter = Array{Float64}(N_bus, N_time, N_scen, N_iter)
a_iter = Array{Float64}(N_bus, N_time, N_scen, N_iter)
pg_iter = Array{Float64}(N_bus, N_time, N_scen, N_iter)
qg_iter = Array{Float64}(N_bus, N_time, N_scen, N_iter)
IC_iter = Array{Float64}(N_scen, N_iter)
OC_iter = Array{Float64}(N_scen, N_iter)
EC_iter = Array{Float64}(N_scen, N_iter)
OCinter_iter = Array{Float64}(N_scen, N_iter)
OCgen_iter = Array{Float64}(N_scen, N_iter)
OCloadshed_iter = Array{Float64}(N_scen, N_iter)
OCstorage_iter = Array{Float64}(N_scen, N_iter)
x_iter = Array{Float64}(N_storage, N_scen, N_iter)
sigma_d_iter = Array{Float64}(N_bus, N_time, N_scen, N_iter)
sigma_l_iter = Array{Float64}(N_line, N_time, N_scen, N_iter)
#########################
penalty1 = Array{Float64}(N_scen,1)
penalty2 = Array{Float64}(N_scen,1)
z_old = Array{Int64}(N_storage, N_bus, N_scen)
x_old = Array{Int64}(N_storage, N_scen)
w_old = PHmultiplier(zeros(N_storage, N_bus, N_scen),zeros(N_storage, N_scen))
w_new = PHmultiplier(zeros(N_storage, N_bus, N_scen),zeros(N_storage, N_scen))

#####################
convergence = Array{Float64}(N_iter,1)
convergence[:] = 777
penalty1_iter = Array{Float64}(N_scen, N_iter)
penalty2_iter = Array{Float64}(N_scen, N_iter)
z_bar_iter = Array{Float64}(N_storage, N_bus, N_iter)
x_bar_iter = Array{Float64}(N_storage, N_iter)
w_iter = PHmultiplier(zeros(N_storage, N_bus, N_scen, N_iter),zeros(N_storage, N_scen, N_iter))

#######
z_old = z_ini
x_old = x_ini
z_bar = sum( Prob[s] * z_old[:,:,s] for s in 1:N_scen)
x_bar = sum( Prob[s] * x_old[:,s] for s in 1:N_scen)


for kk = 1:N_iter
    println("###iteration begin###")
    for s = 1:N_scen
        w_new.z[:,:,s] = w_old.z[:,:,s] + rho.z .* (z_old[:,:,s] - z_bar[:,:] )
        w_new.x[:,s] = w_old.x[:,s] + rho.x .* (x_old[:,s] - x_bar[:] )
        w_new_temp = PHmultiplier(w_new.z[:,:,s],w_new.x[:,s])
        if s in S_N
            contin_temp = contin[:,:,s]
            scen_data_temp = ScenarioData(contin_temp,load24[s,:])
        elseif s in S_E
            contin_temp = contin[:,:,s]
            scen_data_temp = ScenarioData(contin_temp,mean(load24[S_N,:],1))
        end
        status[s],obj_value[s],penalty1[s],penalty2[s],simulation_time[s],ESlocation[:,:,s],SoC[:,:,s],u[:,:,:,s],z[:,:,s],v[:,:,s],a[:,:,s],pg[:,:,s],qg[:,:,s],IC[s],OC[s],EC[s],OCinter[s],OCgen[s],OCloadshed[s],OCstorage[s],x[:,s],sigma_d[:,:,s],sigma_l[:,:,s]  = PHsub_iter(scen_data_temp, gap, delayMode, z_bar, w_new_temp, rho, filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP)
    end
    ###
    z_old = z
    x_old = x
    z_bar = sum( Prob[s] * z_old[:,:,s] for s in 1:N_scen)
    x_bar = sum( Prob[s] * x_old[:,s] for s in 1:N_scen)

    penalty1_iter[:,kk] = penalty1[:]
    penalty2_iter[:,kk] = penalty2[:]
    z_bar_iter[:,:,kk] = z_bar[:,:]
    x_bar_iter[:,kk] = x_bar[:]
    w_iter.z[:,:,:,kk] = w_new.z[:,:,:]
    w_iter.x[:,:,kk] = w_new.x[:,:]


    status_iter[:,kk] = status[:]
    obj_value_iter[:,kk] = obj_value[:]
    simulation_time_iter[:,kk] = simulation_time[:]
    ESlocation_iter[:,:,:,kk] = ESlocation[:,:,:]
    SoC_iter[:,:,:,kk] = SoC[:,:,:]
    u_iter[:,:,:,:,kk] = u[:,:,:,:]
    z_iter[:,:,:,kk] = z[:,:,:]
    v_iter[:,:,:,kk] = v[:,:,:]
    a_iter[:,:,:,kk] = a[:,:,:]
    pg_iter[:,:,:,kk] = pg[:,:,:]
    qg_iter[:,:,:,kk] = qg[:,:,:]
    IC_iter[:,kk] = IC[:]
    OC_iter[:,kk] = OC[:]
    EC_iter[:,kk] = EC[:]
    OCinter_iter[:,kk] = OCinter[:]
    OCgen_iter[:,kk] = OCgen[:]
    OCloadshed_iter[:,kk] = OCloadshed[:]
    OCstorage_iter[:,kk] = OCstorage[:]
    x_iter[:,:,kk] = x[:,:]
    sigma_d_iter[:,:,:,kk] = sigma_d[:,:,:]
    sigma_l_iter[:,:,:,kk] = sigma_l[:,:,:]
    ###

    w_old = w_new
    z_temp = Array{Int64}(N_storage , N_bus, N_scen)
    z_temp[:]=0
    for s = 1:N_scen
        z_temp[:,:,s] = z[:,:,s][:]
    end
    temp_conv = sum( Prob[s] * norm( z_temp[:,:,s][:] - z_bar[:,:][:] ) for s in 1:N_scen)
    obj_value_global_temp = sum(Prob[s]*obj_value[s] for s=1:N_scen)
    println("--------------------")
    println("iteration #:")
    display(kk)
    println()
    println("convergence:")
    display(temp_conv)
    println()
    println("ESlocation:")
    display(ESlocation)
    println()
    println("Objective values:")
    display(obj_value)
    println()
    println("Global objective value:")
    display(obj_value_global_temp)
    println()
    println("Simulation time:")
    display(simulation_time)
    println()
    convergence[kk] = temp_conv
    println("###iteration end###")
    println("--------------------")
end
###########
obj_value_global = Array{Float64}(N_iter,1)
obj_value_global = sum( Prob[s] * obj_value_iter[s,:] for s = 1:N_scen)
simulation_time_series = sum(simulation_time_iter)
simulation_time_parallel=sum(maximum(simulation_time_iter,1))
simulation_time_longest = maximum(simulation_time_iter)

@save "result.jld"
