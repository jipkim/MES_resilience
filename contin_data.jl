function contin_data(N_line,N_time,N_scen_normal,N_scen_contin,N_event)
    # contin: Line Capacity Limiting Coefficient 0~1| normal=0, full trip=1
    contin_normal = Array{Float64}(N_line,N_time,N_scen_normal)
    contin_normal[:] = 0
    if N_event == 1
        contin_7scen = Array{Float64}(N_line,N_time,7)
        contin_7scen[1:end] = 0
        contin_7scen[1,6:end,1] = 0.9
        contin_7scen[12,6:end,2] = 0.9
        contin_7scen[4,6:end,3] = 0.9
        contin_7scen[8,6:end,4] = 0.9
        contin_7scen[9,6:end,5] = 0.9
        contin_7scen[5,6:end,6] = 0.9
        contin_7scen[13,6:end,7] = 0.9

        contin_5scen = Array{Float64}(N_line,N_time,5)
        contin_5scen[1:end] = 0
        contin_5scen[1,6:end,1] = 0.9
        contin_5scen[4,6:end,2] = 0.9
        contin_5scen[5,6:end,3] = 0.9
        contin_5scen[12,6:end,4] = 0.9
        contin_5scen[13,6:end,5] = 0.9

        contin_1scen = Array{Float64}(N_line,N_time,1)
        contin_1scen[1:end] = 0
        contin_1scen[13,6:end,1] = 0.9
    elseif N_event == 2
        contin_7scen = Array{Float64}(N_line,N_time,7)
        contin_7scen[1:end] = 0
        contin_7scen[1,6:end,1] = 0.9
        contin_7scen[12,18:end,1] = 0.9
        contin_7scen[12,6:end,2] = 0.9
        contin_7scen[4,18:end,2] = 0.9
        contin_7scen[4,6:end,3] = 0.9
        contin_7scen[12,18:end,3] = 0.9
        contin_7scen[8,6:end,4] = 0.9
        contin_7scen[12,18:end,4] = 0.9
        contin_7scen[9,6:end,5] = 0.9
        contin_7scen[12,18:end,5] = 0.9
        contin_7scen[5,6:end,6] = 0.9
        contin_7scen[12,18:end,6] = 0.9
        contin_7scen[13,6:end,7] = 0.9
        contin_7scen[4,18:end,7] = 0.9

        contin_5scen = Array{Float64}(N_line,N_time,5)
        contin_5scen[1:end] = 0
        contin_5scen[1,6:end,1] = 0.9
        contin_5scen[8,6:end,1] = 0.9
        contin_5scen[1,6:end,2] = 0.9
        contin_5scen[9,6:end,2] = 0.9
        contin_5scen[4,6:end,3] = 0.9
        contin_5scen[12,6:end,3] = 0.9
        contin_5scen[4,6:end,4] = 0.9
        contin_5scen[13,6:end,4] = 0.9
        contin_5scen[3,6:end,5] = 0.9
        contin_5scen[10,6:end,5] = 0.9

        contin_1scen = Array{Float64}(N_line,N_time,1)
        contin_1scen[1:end] = 0
        contin_1scen[4,6:end,1] = 0.9
        contin_1scen[12,18:end,1] = 0.9

    elseif N_event == 3
        contin_5scen = Array{Float64}(N_line,N_time,5)
        contin_5scen[1:end] = 0

        contin_5scen[1,6:end,1] = 0.9
        contin_5scen[8,6:end,1] = 0.9
        contin_5scen[13,6:end,1] = 0.9

        contin_5scen[1,6:end,2] = 0.9
        contin_5scen[9,6:end,2] = 0.9
        contin_5scen[12,6:end,2] = 0.9

        contin_5scen[1,6:end,3] = 0.9
        contin_5scen[4,6:end,3] = 0.9
        contin_5scen[11,6:end,3] = 0.9

        contin_5scen[4,6:end,4] = 0.9
        contin_5scen[8,6:end,4] = 0.9
        contin_5scen[12,6:end,4] = 0.9

        contin_5scen[1,6:end,5] = 0.9
        contin_5scen[4,6:end,5] = 0.9
        contin_5scen[12,6:end,5] = 0.9
    end



    contin = Array{Float64}(N_line,N_time,N_scen_normal+N_scen_contin)
    contin[:] = 0
    if N_scen_contin == 7
        contin[:,:,N_scen_normal+1:end] = contin_7scen
    elseif N_scen_contin == 5
        contin[:,:,N_scen_normal+1:end] = contin_5scen
    elseif N_scen_contin == 1
        contin[:,:,N_scen_normal+1:end] = contin_1scen
    elseif N_scen_contin == 0
        contin = contin_normal
    else
        error("Thre is no preset contingency scenario data: go edit <contin_data.jl>")
    end

    return contin
end
