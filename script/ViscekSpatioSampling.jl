include("../src/Chemo-MT.jl")


rho = 1
L = 30
j = 1.0
eta = 1.0
N = round(Int, L^2*rho)
# alpha = eta - J*rho
# alpha = -1
# j = 
# j = (eta - alpha) / rho
J = j / π
# # @show alpha

paras = initParas(
                Dr = eta, 
                vel=1.0, 
                J=j, 
                epsilon=0.0, 
                cutoff=2, 
                L=L,
                dx=0.1, 
                sigma=1)

            
spatio_func(x, y) = sin(y * π / L) + 1.0

Jmax = maximum(spatio_func.(0:L, 0:L))*π
Jmin = minimum(spatio_func.(0:L, 0:L))*π

eta_list = 0.0:0.2: (Jmax + 2.0)
N_sample = 20

# save_root = "/mnt/SSD_RAID0/kyle/viscek_spatio/" 
save_root = "Data/"
# save_root = "E:/Kyle/Simulation_Data/"
fname = "VS_sin05_v1/"

save_dir = save_root*fname
function phase_transition_spatio(eta_list, paras, spatio_func; save_dir="")
    #* check if folder exist
    # op_list = similar(eta_list)
    LEN = length(eta_list)
    @show LEN
    for i in 1:LEN
        eta = round(eta_list[i], digits=3)
        @show eta
        paras.Dr = eta
        system = initSystem(paras, N=N)

        @time all_pos, all_ori = simulation_visceks_spatio(system, paras, spatio_func;
            nsteps=50_000,
            dt=0.001,
            isave=60)
        
        if length(save_dir) != 0
            data = Dict("all_pos"=>all_pos, "all_ori"=>all_ori)
            @time save(save_dir*"eta_$(eta).jld2", data)
        end
        # op = get_order_parameter(all_ori)
        # println("ordre parameter = $op")
    end
end

function pt_sampling(N_sample, eta_list; save_dir="")
    # L = 20
    # Jmax = maximum(spatio_func.(0:L, 0:L)) * π
    # Jmin = minimum(spatio_func.(0:L, 0:L)) * π

    # eta_list = 0.0:0.2:(Jmax+2.0)
    p = Progress(N_sample; dt=1.0)
    for Nsamp in 1:N_sample
        save_dir_ = save_dir*"sample$(Nsamp)_"
        @show save_dir_
        @time phase_transition_spatio(eta_list, paras, spatio_func; save_dir=save_dir_)
        # data = Dict("op" => op_list, "eta" => eta_list)
        # save("Data/viscek_spatio/phase_transistion_spatio_05_sample_$i.jld2", data)
        next!(p)
    end
end

pt_sampling(N_sample, eta_list; save_dir=save_dir)


