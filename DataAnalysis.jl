using JLD
# @load "BF-Journal_result_6scen_delayOff.jld"
@load "BF-Journal_result_6scen_delayOn.jld"

status
obj_value
simulation_time
ESlocation
SoC
u
z
pg
qg
pls
qls
IC
OC
EC
OCgen
OCloadshed
OCstorage
x
delta
mu
sigma_d
sigma_l

N_storage = size(u,1)
N_bus = size(u,2)
N_time = size(u,3)
N_scen = size(u,4)


ESlocation_normal = Array{Int64}(N_storage)
temp_z = z
temp_z[find(temp_z.<0.5)]=0
temp_z[find(temp_z.>0.5)]=1
temp_z_idx = ind2sub(temp_z,find(temp_z.==1))
for k = 1:N_storage
	ESlocation_normal[k] = temp_z_idx[2][k] - 1
end

ESlocation = Array{Int64}(N_storage, N_time, N_scen)
for s = 1:N_scen
	for k = 1:N_storage
		temp = u[k,:,:,s]
		temp[find(temp.<0.5)]=0
		temp[find(temp.>0.5)]=1
		temp_idx = ind2sub(temp,find(temp.==1))
		temp_location = Array{Int64}(N_time,1)
		temp_location[1:end] = -1
		for i = 1:length(temp_idx[1])
			temp_location[temp_idx[2][i]]  = temp_idx[1][i] - 1
		end
		ESlocation[k,:,s] = temp_location
	end
end

firstbus = Array{Int64}(N_storage, N_time, N_scen)
lastbus = Array{Int64}(N_storage, N_time, N_scen)
for k in 1:N_storage
	for t in 1:N_time
		for s in 1:N_scen
			firstbus[k,t,s] = findfirst(delta[k,:,t,s],1) - 1
			lastbus[k,t,s] = findlast(delta[k,:,t,s],1) - 1
		end
	end
end

sigma_d[b,t,s];
sigma_l[l,t,s];
ind2sub(sigma_d,find(x->x==0,sigma_d))
ind2sub(sigma_d[:,:,2],find(x->x==0,sigma_d[:,:,2]))
ind2sub(sigma_l,find(x->x==0,sigma_l))







writecsv("result_sigma_d.csv", sigma_d)
