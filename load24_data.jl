## source from: http://www.pjm.com/markets-and-operations/ops-analysis/historical-load-data.aspx
function load24_data(filename_Load24)
    temp_data = readcsv(filename_Load24,header=true)[1][:,3:end]
    raw_data = Array{Float64}(size(temp_data))
    normalized_data = Array{Float64}(size(temp_data))
    raw_data[:] = temp_data[:]
    N_days = size(raw_data,1)
    for k = 1:N_days
        normalized_data[k,:] = raw_data[k,:] / maximum(raw_data[k,:])
    end
    load24 = normalized_data
    return load24
end
