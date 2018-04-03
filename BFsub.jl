function BFsub(scen_data, gap, delayMode, filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP)

    m = Model(solver = GurobiSolver(MIPGap=gap))
    buses, generators, lines, storages, stations, SMP = DataImport(filename_Node, filename_Generator, filename_Line, filename_Storage, filename_Station, filename_SMP)


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
    N_scen = size(scen_data.contin,3)
    S_N = Int64[]
    S_E = Int64[]
    for s = 1:N_scen
        if sum(scen_data.contin,1:2)[s] == 0
            push!(S_N,s)
        elseif sum(scen_data.contin,1:2)[s] > 0
            push!(S_E,s)
        end
    end

    ### 24 hour load profile ###
    # load24 = readcsv("load24.csv")[2,:]
    load24 = Array{Float64}(N_scen,N_time)
    for s in S_N
        load24[s,:] = scen_data.load24[s,:]
    end
    for s in S_E
        load24[s,:] = mean(scen_data.load24,1)
    end


    delay = DelayMatrix(filename_Line)
    t_delay = Array{Int64}(N_bus,N_bus,N_time)
    t_delay_ini = Array{Int64}(N_bus,N_bus)
    for t = 1:N_time
        t_delay[:,:,t] = delay
    end
    t_delay_ini = delay
    #t_delay[1:end] = 2
    #t_delay_ini[1:end] = 2


    b = 1:N_bus
    l = 1:N_line
	i = 1:N_gen
	mg = 1:N_microgrid
    k = 1:N_storage
    h = 1:length(stations)
    t = 1:N_time
    s = 1:N_scen

    PF_ES = tan(acos(0.9))
    PF_LS = tan(acos(0.9))

    Prob = Array{Float64}(N_scen,1) # Scenario Probability
    Prob[S_N] = 0.9 / length(S_N)
    Prob[S_E] = 0.1 / length(S_E)

    LScost = Array{Float64}(1,N_bus)
    LScost[1:end] = 5000
    LScost[13:15] = 10000

    @variable(m, v[b,t,s] >= 0) # variable for voltage square, voltage^2
    @variable(m, a[b,t,s] >= 0) # variable for current

    @variable(m, fp[b,t,s])
    @variable(m, fq[b,t,s])

    @variable(m, pg[b,t,s] >= 0)
    @variable(m, qg[b,t,s])

    @variable(m, e[k,t,s] >= 0)
    @variable(m, SoC[k,t,s] >= 0)
    @variable(m, pch[k,b,t,s] >= 0) # P^ch_(k,b,t): Charging Power of Energy Storage
    @variable(m, pdisch[k,b,t,s] >= 0) # P^disch_(k,b,t): Discharging Power of Energy Storage
    @variable(m, qch[k,b,t,s] >= 0) # P^ch_(k,b,t): Charging Power of Energy Storage
    @variable(m, qdisch[k,b,t,s] >= 0) # P^disch_(k,b,t): Discharging Power of Energy Storage

    @variable(m, u[k,b,t,s], Bin)
    @variable(m, z5[k,t,s], Bin) # ||Journal|| z -> z5 # prevent simantaneous charging
    @variable(m, x[k], Bin)

    @variable(m, z[k,b], Bin) # ||Journal|| Incidence matrix variable of ES unit in stationary mode: 1 if ES unit k is located at bus b, 0 otherwise.

    @variable(m, delta[mg,b,t,s], Bin) # ||Journal||
    @variable(m, mu[mg,l,t,s], Bin) # ||Journal||
    @variable(m, sigma_l[l,t,s], Bin) # ||Journal||
    @variable(m, sigma_d[b,t,s], Bin) # ||Journal||

    @variable(m, IC >= 0)
    @variable(m, OC[s] >= 0)
    @variable(m, EC[s] >= 0) # ||Journal||
    @variable(m, OCgen[s] >= 0)
    @variable(m, OCstorage[s] >= 0)
    @variable(m, OCloadshed[s] >= 0)
	@variable(m, OCinter[s] >= 0)
    @variable(m, PreventMultiSol >= 0)

    #@objective(m, Min, 0)
    @objective(m, Min, IC + sum( Prob[s] * OC[s] for s in S_N) + sum( Prob[s] * EC[s] for s in S_E) + PreventMultiSol)
    #@objective(m, Min, IC + sum( Prob[s] * OC[s] for s in 1:N_scen))
    @constraint(m, PreventMultiSol == 1e-12*sum(sum(sum(sum(b*u[k,b,t,s] for t in 1:N_time) for k in 1:N_storage) for b in 1:N_bus )for s in 1:N_scen))

    DCRF=3.5481*1e-4 # h=10years
    #@constraint(m, IC == sum(DCRF*(storages[k].Cp * storages[k].Pkmax + storages[k].Ce * storages[k].Ekmax) for k=1:N_storage) )
    @constraint(m, IC == sum(DCRF*(storages[k].Cp * storages[k].Pkmax + storages[k].Ce * storages[k].Ekmax) * x[k] for k=1:N_storage) )
    #@constraint(m, IC == sum(storages[k].Cp * storages[k].Pkmax + storages[k].Ce * storages[k].Ekmax for k=1:N_storage) + sum(stations[h].Cf for h=1:length(stations)))
    @constraint(m, OperatingCost[s in S_N], OC[s] == OCgen[s] + OCstorage[s] + OCinter[s])
    @constraint(m, EmergencyCost[s in S_E], EC[s] == OCgen[s] + OCstorage[s] + OCinter[s] + OCloadshed[s] )
    @constraint(m, InterfaceLineCost[s in 1:N_scen], OCinter[s] == sum( SMP[t] * pg[1,t,s] for t = 1:N_time ))
    @constraint(m, GenerationOC[s in 1:N_scen], OCgen[s] == sum(sum(generators[b].cost * pg[b,t,s] for t =1:N_time) for b=2:N_bus))
    @constraint(m, StorageOC[s in 1:N_scen], OCstorage[s] == sum(sum(sum(storages[k].Co * (pch[k,b,t,s] + pdisch[k,b,t,s]) for k=1:N_storage) for b=1:N_bus) for t=1:N_time))
    @constraint(m, LoadShedOC[s in S_E], OCloadshed[s] == sum(sum(LScost[b] * (1-sigma_d[b,t,s]) * buses[b].Pd * load24[s,t] for t =1:N_time) for b=1:N_bus))


    # Energy Storage Initial & Final SoC constraints
    for k = 1:N_storage
        for s in S_N
            @constraint(m, e[k,1,s] == 0.5*storages[k].Ekmax + sum(pch[k,b,1,s] * u[k,b,1,s] * storages[k].Aleph for b=1:N_bus) - sum(pdisch[k,b,1,s] * u[k,b,1,s]/ storages[k].Aleph for b=1:N_bus))
            @constraint(m, e[k,24,s] >= 0.5*storages[k].Ekmax)
        end
        for s in S_E
            @constraint(m, e[k,1,s] == 0.5*storages[k].Ekmax + sum(pch[k,b,1,s] * u[k,b,1,s] * storages[k].Aleph for b=1:N_bus) - sum(pdisch[k,b,1,s] * u[k,b,1,s]/ storages[k].Aleph for b=1:N_bus))
            @constraint(m, e[k,24,s] >= 0*storages[k].Ekmax)
        end
    end

    # Energy Storage Dynamic Equation
    for s in 1:N_scen
        for k = 1:N_storage
            for t = 2:N_time
                @constraint(m, e[k,t,s] == e[k,t-1,s]+sum(pch[k,b,t,s] * u[k,b,t,s] * storages[k].Aleph for b=1:N_bus) - sum(pdisch[k,b,t,s] * u[k,b,t,s]/ storages[k].Aleph for b=1:N_bus))
            end
            for t = 1:N_time
                @constraint(m, SoC[k,t,s] == e[k,t,s] / storages[k].Ekmax)
                @constraint(m, e[k,t,s] <= storages[k].Ekmax)
                for b = 1:N_bus
                    @constraint(m, pch[k,b,t,s]*storages[k].Aleph <= storages[k].Pkmax*z5[k,t,s])
                    @constraint(m, pdisch[k,b,t,s]/storages[k].Aleph <= storages[k].Pkmax*(1-z5[k,t,s]))
                end
            end
        end

        for k = 1:N_storage
            for t = 1:N_time
                for b = 1:N_bus
                    @constraint(m, qch[k,b,t,s] >= - pch[k,b,t,s] * tan(acos(PF_ES)))
                    @constraint(m, qch[k,b,t,s] <= + pch[k,b,t,s] * tan(acos(PF_ES)))
                    @constraint(m, qdisch[k,b,t,s] >= - pdisch[k,b,t,s] * tan(acos(PF_ES)))
                    @constraint(m, qdisch[k,b,t,s] <= + pdisch[k,b,t,s] * tan(acos(PF_ES)))
                end
            end
        end

        for l = 1:N_line
            for t = 1:N_time
                @constraint(m, (fp[lines[l].head,t,s])^2 + (fq[lines[l].head,t,s])^2 <= ( lines[lines[l].arcID].u)^2 )
                @constraint(m, (fp[lines[l].head,t,s] - a[lines[l].head,t,s]*lines[l].r)^2 + (fq[lines[l].head,t,s] - a[lines[l].head,t,s]*lines[l].x)^2 <= ( lines[lines[l].arcID].u)^2 )
                if scen_data.contin[l,t,s] == 0
                    @constraint(m, v[lines[l].head,t,s] - 2*(lines[l].r*fp[lines[l].head,t,s] + lines[l].x*fq[lines[l].head,t,s]) + a[lines[l].head,t,s]*(lines[l].r^2 + lines[l].x^2) == v[lines[l].tail,t,s])
                end
                @constraint(m, ((fp[lines[l].head,t,s])^2 + (fq[lines[l].head,t,s])^2) <= v[lines[l].head,t,s] * a[lines[l].head,t,s] ) # SOCP CONSTRAINT
            end
        end

        for t = 1:N_time
            @constraint(m, v[buses[1].root,t,s] == 1) #Voltage constraint for root node
            @constraint(m, a[buses[1].root,t,s] == 0)
            @constraint(m, fp[buses[1].root,t,s] == 0) # one-way flow: only from Transmission to Distribution
            @constraint(m, fq[buses[1].root,t,s] == 0)
        end
    end

    for s in S_N
        for b = 1:N_bus
            for t = 1:N_time
                @constraint(m, fp[b,t,s]  - sum(fp[l,t,s] - buses[l].R*a[l,t,s] for l in buses[b].children) - pg[b,t,s] - sum(pdisch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + sum(pch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + buses[b].Pd * load24[s,t] == 0)
                @constraint(m, fq[b,t,s]  - sum(fq[l,t,s] - buses[l].X*a[l,t,s] for l in buses[b].children) - qg[b,t,s] - sum(qdisch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + sum(qch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + buses[b].Qd * load24[s,t] - v[b,t,s]*buses[b].B == 0)
            end
        end
    end
    for s in 1:N_scen
        for b = 1:N_bus
            for t = 1:N_time
                if s in S_N
                    @constraint(m, v[b,t,s] <= buses[b].Vmax) #upper limit constraint for voltage square
                    @constraint(m, v[b,t,s] >= buses[b].Vmin) #lower limit constraint for voltage square
                end
            end
        end

        for b = 1:N_bus
            for t = 1:N_time
                @constraint(m, pg[b,t,s] <= generators[b].Pgmax)
                @constraint(m, qg[b,t,s] >= -generators[b].Qgmax)
                @constraint(m, qg[b,t,s] <= generators[b].Qgmax)
            end
        end

        for k = 1:N_storage
            #@constraint(m, sum(u[k,b,1,s] for b in 1:N_bus) == 1)
            #@constraint(m, sum(u[k,b,1,s] for b in 1:N_bus) == x[k])
            for t = 1:N_time
                #@constraint(m, sum(u[k,b,t,s] for b in 1:N_bus) <= 1)
                @constraint(m, sum(u[k,b,t,s] for b in 1:N_bus) <= x[k])
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
                @constraint(m, sum( u[k,b,t,s] for k in 1:N_storage) <= N_station[b] )
            end
        end
        # charging station installed location only
        for k in 1:N_storage
            for b in B_not_ch
                for t in 1:N_time
                    @constraint(m, u[k,b,t,s] == 0)
                end
            end
        end
    end

    for k in 1:N_storage
        for b in 1:N_bus
            for t in 1:N_time
                for s in S_N
                    @constraint(m, u[k,b,t,s] == z[k,b] )
                end
            end
        end
    end
    for k in 1:N_storage
        @constraint(m, sum(z[k,b] for b in 1:N_bus) == x[k])
    end

    ##  Linkage z[k,b] between Normal Operation and Emergency Operations
    for s in S_E
        for k in 1:N_storage
            for i in 1:N_bus
                for j in 1:N_bus
                    if (t_delay_ini[i,j]!=0)
                        for tt in 1:t_delay_ini[i,j]
                            if (i!=j) & (tt <= N_time)
                                @constraint(m, z[k,i] <= 1 - u[k,j,tt,s])
                            end
                            if (i!=j) & (1+tt <= N_time)
                                @constraint(m, z[k,i] - u[k,i,1,s] <= 1 - u[k,j,1+tt,s])
                            end
                        end
                    end
                end
            end
        end
    end

    # Delay Model
    for s in S_E
        for k in 1:N_storage
            for i in 1:N_bus
                for j in 1:N_bus
                    for t in 1:N_time
                        if (t_delay[i,j,t]!=0)
                            for tt in 1:t_delay[i,j,t]
                                if (i!=j) & (t+tt <= N_time) & (delayMode==1)
                                    @constraint(m, u[k,i,t,s] - u[k,i,t+1,s] <= 1 - u[k,j,t+tt,s])
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Microgrid operation mode
    for s in S_E
        for l = 1:N_line
            for t = 1:N_time
                @constraint(m, (fp[lines[l].head,t,s])^2 + (fq[lines[l].head,t,s])^2 <= ( sigma_l[l,t,s] * lines[lines[l].arcID].u)^2 )
                @constraint(m, (fp[lines[l].head,t,s] - a[lines[l].head,t,s]*lines[l].r)^2 + (fq[lines[l].head,t,s] - a[lines[l].head,t,s]*lines[l].x)^2 <= ( sigma_l[l,t,s] * lines[lines[l].arcID].u)^2 )
            end
        end
        for b = 1:N_bus
            for t = 1:N_time
                @constraint(m, fp[b,t,s]  - sum(fp[l,t,s] - buses[l].R*a[l,t,s] for l in buses[b].children) - pg[b,t,s] - sum(pdisch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + sum(pch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + sigma_d[b,t,s] * buses[b].Pd * load24[t] == 0)
                @constraint(m, fq[b,t,s]  - sum(fq[l,t,s] - buses[l].X*a[l,t,s] for l in buses[b].children) - qg[b,t,s] - sum(qdisch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + sum(qch[k,b,t,s] * u[k,b,t,s] for k in 1:N_storage) + sigma_d[b,t,s] * buses[b].Qd * load24[t] - v[b,t,s]*buses[b].B == 0)
            end
        end
    end


    ###### left just in case of PH algorithm ######
    # linkage among scenarios
    for k in 1:N_storage
        for b in 1:N_bus
            for s in 2:N_scen
            #    @constraint(m, u[k,b,1,s] == u[k,b,1,1])
            end
        end
    end
    # equation () | Stationary ESS location for normal operation
    for k in 1:N_storage
        for b in 1:N_bus
            for t in 2:N_time
            #    @constraint(m, u[k,b,t,1] == u[k,b,1,1])
            end
        end
    end
    ##############################################
    for t in 1:N_time
        for s in S_N
            for l in 1:N_line
                @constraint(m, sigma_l[l,t,s] == 1)
            end
            for b in 1:N_bus
                @constraint(m, sigma_d[b,t,s] == 1)
            end
        end
    end

    # Microgrid formation
    for t in 1:N_time
        for s in S_E
			temp_contin = scen_data.contin[:,t,s]
            N_cont = sum(temp_contin.!=0)
            @constraint(m, sum(sigma_l[l,t,s] for l in 1:N_line) >= N_line - max(N_cont,N_storage+N_gen))
            for l in 1:N_line
                if scen_data.contin[l,t,s] != 0
                    @constraint(m, sigma_l[l,t,s] == 0)
                    @constraint(m, fp[lines[l].head,t,s] == 0)
                    @constraint(m, fq[lines[l].head,t,s] == 0)
                    @constraint(m, a[lines[l].head,t,s] == 0)
                end
            end
        end
    end


    tic()
    status = solve(m)
    simulation_time = toc()
    obj_value = getobjectivevalue(m)
    result_v = getvalue(v)[:,:,:] # [b,t,s]
    result_a = getvalue(a)[:,:,:] # [b,t,s]
    result_fp = getvalue(fp)[:,:,:] # [b,t,s]
    result_fq = getvalue(fq)[:,:,:] # [b,t,s]
    result_pg = getvalue(pg)[:,:,:] # [i,t,s]
    result_qg = getvalue(qg)[:,:,:] # [i,t,s]
    result_e = getvalue(e)[:,:,:] # [k,t,s]
    result_SoC = getvalue(SoC)[:,:,:] # [k,t,s]
    result_pch = getvalue(pch)[:,:,:,:] # [k,b,t,s]
    result_pdisch = getvalue(pdisch)[:,:,:,:] # [k,b,t,s]
    result_u = getvalue(u)[:,:,:,:] # [k,b,t,s]

    result_z = getvalue(z)[:,:] # [k,b]
    result_delta = getvalue(delta)[:,:,:,:] # [k,b,t,s]
    result_mu = getvalue(mu)[:,:,:,:] # [k,l,t,s]
    result_sigma_l = getvalue(sigma_l)[:,:,:] # [l,t,s]
    result_sigma_d = getvalue(sigma_d)[:,:,:] # [b,t,s]
    result_EC = getvalue(EC)[:]


    result_IC = getvalue(IC)
    result_OC = getvalue(OC)[:]

    result_OCgen = getvalue(OCgen)[:]
    result_OCloadshed = getvalue(OCloadshed)[:]
    result_OCinter = getvalue(OCinter)[:]
    result_OCstorage = getvalue(OCstorage)[:]
    result_x = getvalue(x)[:] # [k]
    result_ESlocation = Array{Int64}(N_storage, N_time, N_scen)
    for s = 1:N_scen
        for k = 1:N_storage
            temp = result_u[k,:,:,s]
            temp[find(temp.<0.5)]=0
            temp[find(temp.>0.5)]=1
            temp_idx = ind2sub(temp,find(temp.==1))
            temp_location = Array{Int64}(N_time,1)
            temp_location[1:end] = -1
            for i = 1:length(temp_idx[1])
                temp_location[temp_idx[2][i]]  = temp_idx[1][i] - 1
            end
            result_ESlocation[k,:,s] = temp_location
        end
    end

    u = Array{Int64}(N_storage, N_bus, N_time, N_scen)
    u[:] = 0
    for s = 1:N_scen
        for i = 1:N_storage
            for j = 1:N_bus
                for k = 1:N_time
                    if (result_u[i,j,k,s] == 1)
                        u[i,j,k,s] = 1
                    elseif (result_u[i,j,k,s] == 0)
                        u[i,j,k,s] = 0
                    end
                end
            end
        end
    end
    z = Array{Int64}(N_storage, N_bus)
    z[:] = 0
    for k in 1:N_storage
        for b in 1:N_bus
            if (result_z[k,b] == 1)
                z[k,b] = 1
            elseif (result_z[k,b] == 0)
                z[k,b] = 0
            end
        end
    end

    return status, obj_value, simulation_time, result_ESlocation, result_SoC, u, z, result_v, result_a, result_pg, result_qg, result_IC, result_OC, result_EC, result_OCinter, result_OCgen, result_OCloadshed, result_OCstorage, result_x, result_sigma_d, result_sigma_l
end
