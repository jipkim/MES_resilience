function PHsub_ini(scen_data, gap, delayMode, filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP)
    m = Model(solver = GurobiSolver(MIPGap=gap))
    buses, generators, lines, storages, stations, SMP = DataImport(filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP)

    ### 24 hour load profile ###
    load24 = scen_data.load24

    ### Code begins here ###
    N_bus = length(buses)
    N_line = length(lines)
    # N_gen = length(generators)
    N_storage = length(storages)

    B_G = Int64[]
    N_gen = 0
    for b in 1:N_bus
        if generators[b].Pgmax > 0
            N_gen += 1
            push!(B_G,b)
        end
    end
    N_microgrid = N_gen + N_storage

    N_station = Array{Int64}(N_bus,1)
    N_station[1:end] = 0
    for b = 1:N_bus
        for h = 1:length(stations)
            if stations[h].busidx == b
                N_station[b] += 1
            end
        end
    end

    N_time = 24
    if sum(scen_data.contin) == 0
        Contin = 0
    elseif sum(scen_data.contin) > 0
        Contin = 1
    end




    delay = DelayMatrix(filename_Line)
    t_delay = Array{Int64}(N_bus,N_bus,N_time)
    t_delay_ini = Array{Int64}(N_bus,N_bus)
    for t = 1:N_time
        t_delay[:,:,t] = delay
    end
    t_delay_ini = delay

    b = 1:N_bus
    l = 1:N_line
    i = 1:N_gen
    mg = 1:N_microgrid
    k = 1:N_storage
    h = 1:length(stations)
    t = 1:N_time

    PF_ES = tan(acos(0.9))
    PF_LS = tan(acos(0.9))

    LScost = Array{Float64}(1,N_bus)
    LScost[1:end] = 5000
    LScost[13:15] = 10000

    @variable(m, v[b,t] >= 0) # variable for voltage square, voltage^2
    @variable(m, a[b,t] >= 0) # variable for current

    @variable(m, fp[b,t])
    @variable(m, fq[b,t])

    @variable(m, pg[b,t] >= 0)
    @variable(m, qg[b,t])

    @variable(m, e[k,t] >= 0)
    @variable(m, SoC[k,t] >= 0)
    @variable(m, pch[k,b,t] >= 0) # P^ch_(k,b,t): Charging Power of Energy Storage
    @variable(m, pdisch[k,b,t] >= 0) # P^disch_(k,b,t): Discharging Power of Energy Storage
    @variable(m, qch[k,b,t] >= 0) # P^ch_(k,b,t): Charging Power of Energy Storage
    @variable(m, qdisch[k,b,t] >= 0) # P^disch_(k,b,t): Discharging Power of Energy Storage

    @variable(m, u[k,b,t], Bin)
    @variable(m, z5[k,t], Bin) # ||Journal|| z -> z5 # prevent simantaneous charging
    @variable(m, x[k], Bin)

    @variable(m, z[k,b], Bin) # ||Journal|| Incidence matrix variable of ES unit in stationary mode: 1 if ES unit k is located at bus b, 0 otherwise.

    @variable(m, sigma_l[l,t], Bin) # ||Journal||
    @variable(m, sigma_d[b,t], Bin) # ||Journal||

    @variable(m, IC >= 0)
    @variable(m, OC >= 0)
    @variable(m, EC >= 0) # ||Journal||
    @variable(m, OCgen >= 0)
    @variable(m, OCstorage >= 0)
    @variable(m, OCloadshed >= 0)
    @variable(m, OCinter >= 0)
    @variable(m, PreventMultiSol >= 0)
    if Contin == 0
        @objective(m, Min, IC + OC + PreventMultiSol)
    elseif Contin == 1
        @objective(m, Min, IC + EC + PreventMultiSol)
    end

    #@objective(m, Min, IC + sum( Prob[s] * OC[s] for s in 1:N_scen))
    @constraint(m, PreventMultiSol == 1e-12*sum(sum(sum(b*u[k,b,t] for t in 1:N_time) for k in 1:N_storage) for b in 1:N_bus ))

    DCRF=3.5481*1e-4 # h=10years
    #@constraint(m, IC == sum(DCRF*(storages[k].Cp * storages[k].Pkmax + storages[k].Ce * storages[k].Ekmax) for k=1:N_storage) )
    @constraint(m, IC == sum(DCRF*(storages[k].Cp * storages[k].Pkmax + storages[k].Ce * storages[k].Ekmax) * x[k] for k=1:N_storage) )
    #@constraint(m, IC == sum(storages[k].Cp * storages[k].Pkmax + storages[k].Ce * storages[k].Ekmax for k=1:N_storage) + sum(stations[h].Cf for h=1:length(stations)))
    @constraint(m, OperatingCost, OC == OCgen + OCstorage + OCinter)
    @constraint(m, EmergencyCost, EC == OCgen + OCstorage + OCinter + OCloadshed )
    @constraint(m, InterfaceLineCost, OCinter == sum( SMP[t] * pg[1,t] for t = 1:N_time ))
    @constraint(m, GenerationOC, OCgen == sum(sum(generators[b].cost * pg[b,t] for t =1:N_time) for b=2:N_bus))
    @constraint(m, StorageOC, OCstorage == sum(sum(sum(storages[k].Co * (pch[k,b,t] + pdisch[k,b,t]) for k=1:N_storage) for b=1:N_bus) for t=1:N_time))
    @constraint(m, LoadShedOC, OCloadshed == sum(sum(LScost[b] * (1-sigma_d[b,t]) * buses[b].Pd * load24[t] for t =1:N_time) for b=1:N_bus))


    # Energy Storage Initial & Final SoC constraints
    for k = 1:N_storage
        if Contin == 0
            @constraint(m, e[k,1] == 0.5*storages[k].Ekmax + sum(pch[k,b,1] * u[k,b,1] * storages[k].Aleph for b=1:N_bus) - sum(pdisch[k,b,1] * u[k,b,1]/ storages[k].Aleph for b=1:N_bus))
            @constraint(m, e[k,24] >= 0.5*storages[k].Ekmax)
        elseif Contin == 1
            @constraint(m, e[k,1] == 0.5*storages[k].Ekmax + sum(pch[k,b,1] * u[k,b,1] * storages[k].Aleph for b=1:N_bus) - sum(pdisch[k,b,1] * u[k,b,1]/ storages[k].Aleph for b=1:N_bus))
            @constraint(m, e[k,24] >= 0*storages[k].Ekmax)
        end
    end

    # Energy Storage Dynamic Equation

    for k = 1:N_storage
        for t = 2:N_time
            @constraint(m, e[k,t] == e[k,t-1]+sum(pch[k,b,t] * u[k,b,t] * storages[k].Aleph for b=1:N_bus) - sum(pdisch[k,b,t] * u[k,b,t]/ storages[k].Aleph for b=1:N_bus))
        end
        for t = 1:N_time
            @constraint(m, SoC[k,t] == e[k,t] / storages[k].Ekmax)
            @constraint(m, e[k,t] <= storages[k].Ekmax)
            for b = 1:N_bus
                @constraint(m, pch[k,b,t]*storages[k].Aleph <= storages[k].Pkmax*z5[k,t])
                @constraint(m, pdisch[k,b,t]/storages[k].Aleph <= storages[k].Pkmax*(1-z5[k,t]))
            end
        end
    end

    for k = 1:N_storage
        for t = 1:N_time
            for b = 1:N_bus
                @constraint(m, qch[k,b,t] >= - pch[k,b,t] * tan(acos(PF_ES)))
                @constraint(m, qch[k,b,t] <= + pch[k,b,t] * tan(acos(PF_ES)))
                @constraint(m, qdisch[k,b,t] >= - pdisch[k,b,t] * tan(acos(PF_ES)))
                @constraint(m, qdisch[k,b,t] <= + pdisch[k,b,t] * tan(acos(PF_ES)))
            end
        end
    end

    for l = 1:N_line
        for t = 1:N_time
            @constraint(m, (fp[lines[l].head,t])^2 + (fq[lines[l].head,t])^2 <= ( lines[lines[l].arcID].u)^2 )
            @constraint(m, (fp[lines[l].head,t] - a[lines[l].head,t]*lines[l].r)^2 + (fq[lines[l].head,t] - a[lines[l].head,t]*lines[l].x)^2 <= ( lines[lines[l].arcID].u)^2 )
            if Contin == 0
                @constraint(m, v[lines[l].head,t] - 2*(lines[l].r*fp[lines[l].head,t] + lines[l].x*fq[lines[l].head,t]) + a[lines[l].head,t]*(lines[l].r^2 + lines[l].x^2) == v[lines[l].tail,t])
            elseif Contin == 1
                if scen_data.contin[l,t] == 0
                    @constraint(m, v[lines[l].head,t] - 2*(lines[l].r*fp[lines[l].head,t] + lines[l].x*fq[lines[l].head,t]) + a[lines[l].head,t]*(lines[l].r^2 + lines[l].x^2) == v[lines[l].tail,t])
                end
            end
            @constraint(m, ((fp[lines[l].head,t])^2 + (fq[lines[l].head,t])^2) <= v[lines[l].head,t] * a[lines[l].head,t] ) # SOCP CONSTRAINT
        end
    end

    for t = 1:N_time
        @constraint(m, v[buses[1].root,t] == 1) #Voltage constraint for root node
        @constraint(m, a[buses[1].root,t] == 0)
        @constraint(m, fp[buses[1].root,t] == 0) # one-way flow: only from Transmission to Distribution
        @constraint(m, fq[buses[1].root,t] == 0)
    end


    if Contin == 0
        for b = 1:N_bus
            for t = 1:N_time
                @constraint(m, fp[b,t]  - sum(fp[l,t] - buses[l].R*a[l,t] for l in buses[b].children) - pg[b,t] - sum(pdisch[k,b,t] * u[k,b,t] for k in 1:N_storage) + sum(pch[k,b,t] * u[k,b,t] for k in 1:N_storage) + buses[b].Pd * load24[t] == 0)
                @constraint(m, fq[b,t]  - sum(fq[l,t] - buses[l].X*a[l,t] for l in buses[b].children) - qg[b,t] - sum(qdisch[k,b,t] * u[k,b,t] for k in 1:N_storage) + sum(qch[k,b,t] * u[k,b,t] for k in 1:N_storage) + buses[b].Qd * load24[t] - v[b,t]*buses[b].B == 0)
            end
        end
    end

    for b = 1:N_bus
        for t = 1:N_time
            if Contin == 0
                @constraint(m, v[b,t] <= buses[b].Vmax) #upper limit constraint for voltage square
                @constraint(m, v[b,t] >= buses[b].Vmin) #lower limit constraint for voltage square
            end
        end
    end

    for b = 1:N_bus
        for t = 1:N_time
            @constraint(m, pg[b,t] <= generators[b].Pgmax)
            @constraint(m, qg[b,t] >= -generators[b].Qgmax)
            @constraint(m, qg[b,t] <= generators[b].Qgmax)
        end
    end

    for k = 1:N_storage
        #@constraint(m, sum(u[k,b,1] for b in 1:N_bus) == x[k])
        for t = 1:N_time
            @constraint(m, sum(u[k,b,t] for b in 1:N_bus) <= x[k])
        end
    end

    #**** index generation
    temp = collect(1:N_bus)
    for h = 1:length(stations)
        temp_idx = find(temp .== stations[h].busidx)
        if (!isempty(temp_idx))
            splice!(temp,temp_idx[1])
        end
    end
    B_not_ch = temp
    temp = collect(1:N_bus)
    B_ch = setdiff(temp,B_not_ch)

    # equation (19)
    for t in 1:N_time
        for b in B_ch
            @constraint(m, sum( u[k,b,t] for k in 1:N_storage) <= N_station[b] )
        end
    end
    # charging station installed location only
    for k in 1:N_storage
        for b in B_not_ch
            for t in 1:N_time
                @constraint(m, u[k,b,t] == 0)
            end
        end
    end


    if Contin == 0
        for k in 1:N_storage
            for b in 1:N_bus
                for t in 1:N_time
                    @constraint(m, u[k,b,t] == z[k,b] )
                end
            end
        end
    end
    for k in 1:N_storage
        @constraint(m, sum(z[k,b] for b in 1:N_bus) == x[k])
    end


    ##  Linkage z[k,b] between Normal Operation and Emergency Operations
    if Contin == 1
        for k in 1:N_storage
            for i in 1:N_bus
                for j in 1:N_bus
                    if (t_delay_ini[i,j]!=0)
                        for tt in 1:t_delay_ini[i,j]
                            if (i!=j) & (tt <= N_time)
                                @constraint(m, z[k,i] <= 1 - u[k,j,tt])
                            end
                            if (i!=j) & (1+tt <= N_time)
                                @constraint(m, z[k,i] - u[k,i,1] <= 1 - u[k,j,1+tt])
                            end
                        end
                    end
                end
            end
        end
    end

    # Delay Model
    if Contin == 1
        for k in 1:N_storage
            for i in 1:N_bus
                for j in 1:N_bus
                    for t in 1:N_time
                        if (t_delay[i,j,t]!=0)
                            for tt in 1:t_delay[i,j,t]
                                if (i!=j) & (t+tt <= N_time) & (delayMode==1)
                                    @constraint(m, u[k,i,t] - u[k,i,t+1] <= 1 - u[k,j,t+tt])
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Microgrid operation mode
    if Contin == 1
        for l = 1:N_line
            for t = 1:N_time
                @constraint(m, (fp[lines[l].head,t])^2 + (fq[lines[l].head,t])^2 <= ( sigma_l[l,t] * lines[lines[l].arcID].u)^2 )
                @constraint(m, (fp[lines[l].head,t] - a[lines[l].head,t]*lines[l].r)^2 + (fq[lines[l].head,t] - a[lines[l].head,t]*lines[l].x)^2 <= ( sigma_l[l,t] * lines[lines[l].arcID].u)^2 )
            end
        end
        for b = 1:N_bus
            for t = 1:N_time
                @constraint(m, fp[b,t]  - sum(fp[l,t] - buses[l].R*a[l,t] for l in buses[b].children) - pg[b,t] - sum(pdisch[k,b,t] * u[k,b,t] for k in 1:N_storage) + sum(pch[k,b,t] * u[k,b,t] for k in 1:N_storage) + sigma_d[b,t] * buses[b].Pd * load24[t] == 0)
                @constraint(m, fq[b,t]  - sum(fq[l,t] - buses[l].X*a[l,t] for l in buses[b].children) - qg[b,t] - sum(qdisch[k,b,t] * u[k,b,t] for k in 1:N_storage) + sum(qch[k,b,t] * u[k,b,t] for k in 1:N_storage) + sigma_d[b,t] * buses[b].Qd * load24[t] - v[b,t]*buses[b].B == 0)
            end
        end
    end

    if Contin == 0
        for t in 1:N_time
            for l in 1:N_line
                @constraint(m, sigma_l[l,t] == 1)
            end
            for b in 1:N_bus
                @constraint(m, sigma_d[b,t] == 1)
            end
        end
    end

    # Microgrid formation
    if Contin == 1
        for t in 1:N_time
            temp_contin = scen_data.contin[:,t]
            N_cont = sum(temp_contin.!=0)
            @constraint(m, sum(sigma_l[l,t] for l in 1:N_line) >= N_line - max(N_cont,N_storage+N_gen))
            for l in 1:N_line
                if scen_data.contin[l,t] != 0
                    @constraint(m, sigma_l[l,t] == 0)
                    @constraint(m, fp[lines[l].head,t] == 0)
                    @constraint(m, fq[lines[l].head,t] == 0)
                    @constraint(m, a[lines[l].head,t] == 0)
                end
            end
        end
    end


    tic()
    status = solve(m)
    simulation_time = toc()
    obj_value = getobjectivevalue(m)
    result_v = getvalue(v)[:,:] # [b,t]
    result_a = getvalue(a)[:,:] # [b,t]
    result_fp = getvalue(fp)[:,:] # [b,t]
    result_fq = getvalue(fq)[:,:] # [b,t]
    result_pg = getvalue(pg)[:,:] # [i,t]
    result_qg = getvalue(qg)[:,:] # [i,t]
    result_e = getvalue(e)[:,:] # [k,t]
    result_SoC = getvalue(SoC)[:,:] # [k,t]
    result_pch = getvalue(pch)[:,:,:] # [k,b,t]
    result_pdisch = getvalue(pdisch)[:,:,:] # [k,b,t]
    result_u = getvalue(u)[:,:,:] # [k,b,t]

    result_z = getvalue(z)[:,:] # [k,b]
    result_sigma_l = getvalue(sigma_l)[:,:] # [l,t]
    result_sigma_d = getvalue(sigma_d)[:,:] # [b,t]
    result_EC = getvalue(EC)
    result_IC = getvalue(IC)
    result_OC = getvalue(OC)
    result_OCgen = getvalue(OCgen)
    result_OCloadshed = getvalue(OCloadshed)
    result_OCinter = getvalue(OCinter)
    result_OCstorage = getvalue(OCstorage)
    result_x = getvalue(x)[:] # [k]
    result_ESlocation = Array{Int64}(N_storage, N_time)
    for k = 1:N_storage
        temp = result_u[k,:,:]
        temp[find(temp.<0.5)]=0
        temp[find(temp.>0.5)]=1
        temp_idx = ind2sub(temp,find(temp.==1))
        temp_location = Array{Int64}(N_time,1)
        temp_location[1:end] = -1
        for i = 1:length(temp_idx[1])
            temp_location[temp_idx[2][i]]  = temp_idx[1][i] - 1
        end
        result_ESlocation[k,:] = temp_location
    end


    u = Array{Int64}(N_storage, N_bus, N_time)
    u[:] = -1
    for i = 1:N_storage
        for j = 1:N_bus
            for k = 1:N_time
                if (result_u[i,j,k] >= 0.9)
                    u[i,j,k] = 1
                elseif (result_u[i,j,k] <= 0.1)
                    u[i,j,k] = 0
                end
            end
        end
    end

    z = Array{Int64}(N_storage, N_bus)
    z[:] = -1
    for k in 1:N_storage
        for b in 1:N_bus
            if (result_z[k,b] >= 0.9)
                z[k,b] = 1
            elseif (result_z[k,b] <= 0.1)
                z[k,b] = 0
            end
        end
    end
#z,
    return status, obj_value, simulation_time, result_ESlocation, result_SoC, u, z, result_v, result_a, result_pg, result_qg, result_IC, result_OC, result_EC, result_OCinter, result_OCgen, result_OCloadshed, result_OCstorage, result_x, result_sigma_d, result_sigma_l
end
